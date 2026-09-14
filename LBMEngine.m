%{
    LBM
%}

classdef LBMEngine < handle

    properties %(Access=private)

        redblue = [1 0 0; 1 1 1; 0 0 1]

        Nx double; % wind tunnel length
        Ny double; % wind tunnel height
        u0x double = 0.05; % inlet velocity (x only) (constant)
        c = sqrt(1/3);

        %(Ny,Nx)
        rho (:,:) double; % density
        % u (:,:,2) double; % velocities
        ux (:,:) double; % x velocity
        uy (:,:) double; % y velocity

        f (:,:,9) double % particle PDFs
        feq (:,:,9) double;  % equilibium PDFs
        fstar (:,:,9) double; % particle PDFs after collision before stream

        w (1,9) double = [ 4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36 ]; %weights
        % c (1,9,2) double; % direction vectors
        cx (1,9) double = [ 0,  1,  0, -1,  0,  1, -1, -1,  1 ]; % direction vector x
        cy (1,9) double = [ 0,  0,  1,  0, -1,  1,  1, -1, -1 ]; % direction vector y      
        opp (1,9) double = [1, 4, 5, 2, 3, 8, 9, 6, 7]; 

            % 7 3 6
            % 4 1 2
            % 8 5 9

        tau double = 0.55; % relaxation time

        airfoil NACAAirfoil;
        airfoilOrigin (2,1) double; % [x;y]

        obstacle (:,:) logical; % thing to be tested in wind tunnel
        obs3d (:,:,9) logical;

    end

    methods (Access=public)
        function obj = LBMEngine(Nx, Ny)
            arguments
                Nx double; % wind tunnel length
                Ny double; % wind tunnel height
            end

            obj.Nx = Nx;
            obj.Ny = Ny;

            % obj.f = zeros(Ny,Nx,9);
            % obj.f(:,:,1) = ones(Ny,Nx,1);
            % obj.feq = zeros(Ny,Nx,9);

            obj.rho = ones(Ny, Nx);
            obj.ux = ones(Ny, Nx) * obj.u0x;
            % obj.uy = zeros(Ny, Nx);
            obj.uy = repmat(1e-4*sin(2*pi*(1:obj.Ny)'/obj.Ny), 1, Nx);
            obj.f = zeros(Ny, Nx, 9);
            obj.feq = zeros(Ny,Nx,9);
            obj.calcFeq();
            obj.f = obj.feq;



        end

        function addObstacle(obj, obs)
            arguments
                obj 
                obs (:,:) logical;
            end
            obj.obstacle = obs; %repmat(obs, [1 1 9]);
            obj.obs3d = repmat(obs, [1 1 9]);
        end
        function addAirfoil(obj, airfoil, x, y, L, aoa)
            arguments
                obj 
                airfoil NACAAirfoil
                x double
                y double
                L double
                aoa double
            end
            obj.airfoil = airfoil;
            obj.airfoilOrigin = [x;y];
            obj.obstacle = airfoil.makeObstacle(x,y, obj.Nx,obj.Ny, L, aoa);
            obj.obs3d = repmat(obj.obstacle, [1 1 9]);
            obj.Re = obj.u0x*L/obj.nu;
        end


        function runDensity(obj, t)
            arguments
                obj 
                t double; % total time to run (unitless)
            end

            wtp = obj.setupPlot();
            clim(wtp.Parent, [0.95 1.05]);
            
            for T = 0:t
                obj.calcRhoVel();
                obj.calcFeq();
                obj.collision();
                obj.stream();
                obj.applyBC();

                wtp.CData = obj.applyObs(obj.rho, obj.obstacle);
                drawnow limitrate;

                if mod(T,500) == 0
                    disp(T + "/" + t)
                end
            end
            
        end

        function runPressure(obj, t)
            arguments
                obj 
                t double; % total time to run (unitless)
            end

            wtp = obj.setupPlot();
            clim(wtp.Parent, [0.328 0.332]);
            
            for T = 0:t
                obj.calcRhoVel();
                obj.calcFeq();
                obj.collision();
                obj.stream();
                obj.applyBC();

                wtp.CData = obj.applyObs(obj.c^2*obj.rho, obj.obstacle);
                drawnow limitrate;

                if mod(T,500) == 0
                    disp(T + "/" + t)
                end
            end
        end
        
        function runVelocity(obj, t)
            arguments
                obj 
                t double; % total time to run (unitless)
            end

            wtp = obj.setupPlot();
            clim(wtp.Parent, [0.03 0.07]);
            
            for T = 0:t
                obj.calcRhoVel();
                obj.calcFeq();
                obj.collision();
                obj.stream();
                obj.applyBC();

                wtp.CData = obj.applyObs(sqrt(obj.ux.^2 + obj.uy.^2), obj.obstacle);
                drawnow limitrate;

                if mod(T,500) == 0
                    disp(T + "/" + t)
                end
            end
            
        end

        function runVorticity(obj, t)
            arguments
                obj 
                t double; % total time to run (unitless)
            end

            function vort = calcVorticity(ux,uy)
                vort = zeros(obj.Ny,obj.Nx);
                vort(2:end-1, 2:end-1) = 0.5 * (uy(2:end-1, 3:end) - uy(2:end-1, 1:end-2)) ...
                    - 0.5 * (ux(3:end, 2:end-1) - ux(1:end-2, 2:end-1));
                % for i = 2:obj.Ny-1
                %     for j = 2:obj.Nx-1
                %         vort(i,j) = 0.5*(uy(i,j+1)-uy(i,j-1)) - 0.5*(ux(i+1,j)-ux(i-1,j));
                %     end
                % end
            end

            wtp = obj.setupPlot()
            clim(wtp.Parent, [-2e-3 2e-3]);
            
            
            for T = 0:t
                obj.calcRhoVel();
                obj.calcFeq();
                obj.collision();
                obj.stream();
                obj.applyBC();

                % h.CData = applyObs(calcVorticity(obj.ux, obj.uy), obj.obstacle);
                % drawnow limitrate;

                if mod(T,500) == 0
                    disp(T + "/" + t)
                    wtp.CData = obj.applyObs(calcVorticity(obj.ux, obj.uy), obj.obstacle);
                    drawnow;
                end
            end
            
        end
    end

    methods (Access=private)

        function calcRhoVel(obj)
            obj.rho = sum(obj.f,3);
            obj.ux = sum(obj.f .* reshape(obj.cx, 1,1,9),3) ./ obj.rho;
            obj.uy = sum(obj.f .* reshape(obj.cy, 1,1,9),3) ./ obj.rho;
        end

        function calcFeq(obj)
            u2 = obj.ux.^2+obj.uy.^2;
            fEqLocal = obj.feq;

            for i = 1:9
                cdotu = obj.cx(i)*obj.ux+obj.cy(i)*obj.uy;
                fEqLocal(:,:,i) = obj.rho .* obj.w(i) .* ...
                    (1 + (cdotu)/(obj.c^2) ...
                       + (cdotu).^2/(2*obj.c^4) ...
                       - (u2)/(2*obj.c^2));
            end

            obj.feq = fEqLocal;
        end

        function collision(obj)
            obj.fstar = obj.f - (1/obj.tau) * (obj.f-obj.feq);

            bounceBack = obj.f(:,:,obj.opp);
            obs3d = repmat(obj.obstacle, [1 1 9]);
            obj.fstar(obs3d) = bounceBack(obs3d);
        end

        function stream(obj)
            fLocal = obj.f;
            fStarLocal = obj.fstar;
            
            for i = 1:9
                fLocal(:,:,i) = circshift(fStarLocal(:,:,i), [obj.cy(i), obj.cx(i)]);
            end
            obj.f = fLocal;
        end

        function applyBC(obj) % Zou-He
            % inlet (main)
            rho_inlet = (1/(1-obj.u0x)) * (obj.f(:,1,1)+obj.f(:,1,3)+obj.f(:,1,5) + 2*(obj.f(:,1,4)+obj.f(:,1,7)+obj.f(:,1,8)));
            obj.f(:,1,2) = obj.f(:,1,4) + (2/3)*rho_inlet.*obj.u0x;
            obj.f(:,1,6) = obj.f(:,1,8) + (1/6)*rho_inlet.*obj.u0x - (1/2)*(obj.f(:,1,3)-obj.f(:,1,5));
            obj.f(:,1,9) = obj.f(:,1,7) + (1/6)*rho_inlet.*obj.u0x + (1/2)*(obj.f(:,1,3)-obj.f(:,1,5));
            % % inlet (top)
            % obj.f(1,1,2) = obj.f(1,1,4);
            % obj.f(1,1,5) = obj.f(1,1,3);
            % obj.f(1,1,9) = obj.f(1,1,7);
            % obj.f(1,1,6) = (0.5)*(rho_inlet(1) - sum(obj.f(1,1,[1,2,3,4,5,7,9])));
            % obj.f(1,1,8) = obj.f(1,1,6);
            % % inlet (bottom)
            % obj.f(end,1,2) = obj.f(end,1,4);
            % obj.f(end,1,3) = obj.f(end,1,5);
            % obj.f(end,1,6) = obj.f(end,1,8);
            % obj.f(end,1,7) = (0.5)*(rho_inlet(end) - sum(obj.f(end,1,[1,2,3,4,5,6,8])));
            % obj.f(end,1,9) = obj.f(end,1,7);

            % outlet
            % ux_outlet = -1 + obj.f(:,end,1)+obj.f(:,end,3)+obj.f(:,end,5) + 2*(obj.f(:,end,2)+obj.f(:,end,6)+obj.f(:,end,9));
            % obj.f(:,end,4) = obj.f(:,end,2) - (2/3)*ux_outlet;
            % obj.f(:,end,7) = obj.f(:,end,9) - (1/6)*ux_outlet - (1/2)*(obj.f(:,end,3)-obj.f(:,end,5));
            % obj.f(:,end,8) = obj.f(:,end,6) - (1/6)*ux_outlet + (1/2)*(obj.f(:,end,3)-obj.f(:,end,5));
            obj.f(:,end,:)=obj.f(:,end-1,:);


            % walls
            fLocal = obj.f;
            fStarLocal = obj.fstar;
            for i = 1:9
                fLocal(1,:, i) = fStarLocal(1,:, obj.opp(i));
                fLocal(end,:, i) = fStarLocal(end,:, obj.opp(i));
            end
            obj.f = fLocal;

            % obstacle          
            % bounceBack = zeros(obj.Ny,obj.Nx,9);
            % bounceBack(:,:,1:9) = obj.fstar(:,:,obj.opp);
            % obs3d = repmat(obj.obstacle, [1 1 9]);
            % obj.f(obs3d) = bounceBack(obs3d);
            

        end

        
        function wtp = setupPlot(obj)
            figure;
            wtp = imagesc(obj.applyObs(obj.rho, obj.obstacle));
            xlim([0 obj.Nx]); ylim([0 obj.Ny]);
            axis equal;
            colormap("nebula"); colorbar; 
            % clim([0.95 1.05]);
            set(wtp, 'AlphaData', ~obj.obstacle)
        end
    
        function maskedData = applyObs(obj, data, obs)
            maskedData = data;
            maskedData(obs(:,:,1)) = NaN;
        end
    end

    methods (Static)
        function main()
            Ny = 500;
            Nx = 1500;
           
            % obs = nan(Ny,Nx);
            % for i = 1:Ny
            %     for j = 1:Nx
            %         obs(i,j) = (i-Ny*.48)^2 + (j-Nx*.15)^2 <= (Ny*.075)^2;
            %     end
            % end
            % 
            % test = LBMEngine(Nx,Ny);
            % test.addObstacle(obs);
            % test.runVorticity(25000);
            
            test = LBMEngine(Nx, Ny);

            naca4412 = NACAAirfoilObstacle(4412);
            disp("naca aifoil")
            obs = naca4412.makeObstacle(225,225, Nx,Ny, 250, 25);
            disp("makeObstacle")
            test.addObstacle(obs);
            
            disp("run")
            test.runVorticity(25000);

        end


    end

end
