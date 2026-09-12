
classdef NACAAirfoil < handle

    properties (Access=public)
        name double;
        m double;
        p double
        tt double;
        
        L double=1; % length of airfoil
        aoa double=0; % angle of attack
        CL function_handle; % Coefficient of Lift
        
        Axu function_handle;
        Ayu function_handle;
        Axl function_handle;
        Ayl function_handle;
        Au  function_handle;
        Al  function_handle;
        dAxu function_handle;
        dAyu function_handle;
        dAxl function_handle;
        dAyl function_handle;
        dAu  function_handle;
        dAl  function_handle;
    end
    properties (Access=private)
        wu function_handle;
        wl function_handle;
    end

    methods %(Static)
        function obj = NACAAirfoil(name)
            arguments
                name double %{mustBeValid(name)};
            end

            obj.m = floor(name/1000)/100;
            obj.p = mod(floor(name/100),10)/10;
            obj.tt = mod(name,100)/100;

            obj.calculateFunctions();

            if (obj.m == 0 && obj.p == 0)
                obj.p=0.01;
            end
            
        end


        function obs = makeObstacle(obj, x, y, l, h, L, aoa)
            arguments
                obj;
                % (x,y) = (0,0) is top left corner
                x double; % x location of airfoil
                y double; % y location of airfoil
                l double; % length (x) of wind tunnel
                h double; % height (h) of wind tunnel
                L double; % length of airfoil
                aoa double = 0; % angle of attack (degrees)
            end

            obj.aoa = aoa;
            obj.L = L;
            % obj.calculateFunctions();            
            
            rL = ceil(L*cosd(aoa)); % relative length (to the x-axis)
            [X,Y] = meshgrid(x-1:x+rL+1, 1:h);
            X = (X-x); Y = -(Y-y);
            obsu = sign(obj.wu(X,Y));
            obsl = sign(obj.wl(X,Y));
            obs = nan(h,l);
            obs(:,x-1:x+rL+1) = abs(.5*(obsu+obsl)) == 0;
            obs(:,1:x-1) = false;
            obs(:,x+rL+1:l) = false;
            

            % figure;
            % imagesc(obs); colorbar;
            % axis equal;
        end
    end

    methods (Access=public)
        function calculateFunctions(obj)

            yt = @(x) 5*obj.tt*(0.2969*x.^.5 - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);
            dyt = @(x) 5*obj.tt*(0.14845*x.^-.5 - 0.1260 - 0.7032*x + 0.8529*x.^2 - 0.4060*x.^3);
            yc = @(x) ((0<=x&x<=obj.p) .* (obj.m*(2*obj.p*x-x.^2)/obj.p^2)) + ((obj.p<=x & x<=1) .* (obj.m*(1-2*obj.p+2*obj.p*x-x.^2)/(1-obj.p)^2));
            dyc = @(x) ((0<=x&x<=obj.p) .* (2*obj.m*(obj.p-x)/obj.p^2)) + ((obj.p<=x & x<=1) .* (2*obj.m*(obj.p-x)/(1-obj.p)^2));
            d2yc = @(x) ((0<=x&x<=obj.p) .* (-2*obj.m/obj.p^2)) + ((obj.p<=x & x<=1) .* (-2*obj.m/(1-obj.p)^2));

            % aoa = 0
            % xu = @(t) L*(t - yt(t).*dyc(t)./sqrt(1+dyc(t).^2));
            % dxu = @(t) L*(1 - (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5));
            % xl = @(t) L*(t + yt(t).*dyc(t)/sqrt(1+dyc(t).^2));
            % dxl = @(t) L*(1 + (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5));
            % yu = @(t) L*(yc(t) + yt(t)./(1+dyc(t).^2).^.5);
            % dyu = @(t) L*(dyc(t) + (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5));
            % yl = @(t) L*(yc(t) - yt(t)./(1+dyc(t).^2).^.5);
            % dyl = @(t) L*(dyc(t) - (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5));

            xu = @(t) obj.L*(t - yt(t).*dyc(t)./sqrt(1+dyc(t).^2))*cosd(obj.aoa) + obj.L*(yc(t) + yt(t)./(1+dyc(t).^2).^.5)*sind(obj.aoa);
            dxu = @(t) obj.L*(1 - (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*cosd(obj.aoa) ...
                       + obj.L*(dyc(t) + (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*sind(obj.aoa);
            yu = @(t) obj.L*(yc(t) + yt(t)./(1+dyc(t).^2).^.5)*cosd(obj.aoa) - obj.L*(t - yt(t).*dyc(t)./sqrt(1+dyc(t).^2))*sind(obj.aoa);
            dyu = @(t) obj.L*(dyc(t) + (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*cosd(obj.aoa) ...
                       - obj.L*(1 - (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*sind(obj.aoa);
            xl = @(t) obj.L*(t + yt(t).*dyc(t)/sqrt(1+dyc(t).^2))*cosd(obj.aoa) + obj.L*(yc(t) - yt(t)./(1+dyc(t).^2).^.5)*sind(obj.aoa);
            dxl = @(t) obj.L*(1 + (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*cosd(obj.aoa) ...
                       + obj.L*(dyc(t) - (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*sind(obj.aoa);
            yl = @(t) obj.L*(yc(t) - yt(t)./(1+dyc(t).^2).^.5)*cosd(obj.aoa) - obj.L*(t + yt(t).*dyc(t)/sqrt(1+dyc(t).^2))*sind(obj.aoa);
            dyl = @(t) obj.L*(dyc(t) - (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*cosd(obj.aoa) ...
                       - obj.L*(1 + (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*sind(obj.aoa);

            % x, y position of airfoil (parametric wrt t)
            obj.Axu = @(t) xu(t).*cosd(obj.aoa) + yu(t).*sind(obj.aoa);
            obj.dAxu = @(t) dxu(t).*cosd(obj.aoa) + dyu(t).*sind(obj.aoa);
            obj.Ayu = @(t) yu(t).*cosd(obj.aoa) - xu(t).*sind(obj.aoa);
            obj.dAyu = @(t) dyu(t).*cosd(obj.aoa) - dxu(t).*sind(obj.aoa);
            obj.Au = @(t) [obj.Axu(t); obj.Ayu(t)];
            obj.dAu = @(t) [obj.dAxu(t); obj.dAyu(t)];

            obj.Axl = @(t) xl(t).*cosd(obj.aoa) + yl(t).*sind(obj.aoa);
            obj.dAxl = @(t) dxl(t).*cosd(obj.aoa) + dyl(t).*sind(obj.aoa);
            obj.Ayl = @(t) yl(t).*cosd(obj.aoa) - xl(t).*sind(obj.aoa);
            obj.dAyl = @(t) dyl(t).*cosd(obj.aoa) - dxl(t).*sind(obj.aoa);
            obj.Al = @(t) [obj.Axl(t); obj.Ayl(t)];
            obj.dAl = @(t) [obj.dAxl(t); obj.dAyl(t)];

            % find function name
            obj.wu = @(x,y) integral(@(t) ((xu(t)-x).*dyu(t) - (yu(t)-y).*dxu(t))./((xu(t)-x).^2 + (yu(t)-y).^2), 0,1, 'ArrayValued',true);
            obj.wl = @(x,y) integral(@(t) ((xl(t)-x).*dyl(t) - (yl(t)-y).*dxl(t))./((xl(t)-x).^2 + (yl(t)-y).^2), 0,1, 'ArrayValued',true);


            % https://aviation.stackexchange.com/questions/96119/is-there-a-formula-for-calculating-lift-coefficient-based-on-the-naca-airfoil
            thetaP = acos(1-2*obj.p);
            A0 = obj.m*((2*obj.p-1)*thetaP + sin(thetaP))/(pi*obj.p^2) +  obj.m*((2*obj.p-1)*(pi-thetaP) - sin(thetaP))/(pi*(1-obj.p)^2);
            if isnan(A0), A0=0; end
            A1 = 2*obj.m*((2*obj.p-1)*sin(thetaP) + 0.25*sin(2*thetaP) + 0.5*thetaP)/(pi*obj.p^2) - 2*obj.m*((2*obj.p-1)*sin(thetaP) + 0.25*sin(2*thetaP) - 0.5*(pi-thetaP))/(pi*(1-obj.p)^2);
            if isnan(A1), A1=0; end

            obj.CL = @() 2*pi*(obj.aoa*pi/180) + pi*(A1-2*A0);

        end

    %     function mustBeValid(name)
    %         if ~(name>=1000 && name <= 9999)
    %             error("invalid airfoil name. must be 4 digits")
    %         end
    %     end
    end
end