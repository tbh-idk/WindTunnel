
classdef NACAAirfoilObstacle < handle

    properties
        name double;
        m double;
        p double
        tt double;

        wu function_handle;
        wl function_handle;

    end

    methods %(Static)
        function obj = NACAAirfoilObstacle(name)
            arguments
                name double %{mustBeValid(name)};
            end

            obj.m = floor(name/1000)/100;
            obj.p = mod(floor(name/100),10)/10;
            obj.tt = mod(name,100)/100;

            % obj.calculateFunctions(0);
            
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

            obj.calculateFunctions(aoa, L);            

            [X,Y] = meshgrid(x-1:x+L*cosd(aoa)+1, 1:h);
            X = (X-x); Y = -(Y-y);
            obsu = sign(obj.wu(X,Y));
            obsl = sign(obj.wl(X,Y));
            obs = nan(h,l);
            obs(:,x-1:x+L*cosd(aoa)+1) = abs(.5*(obsu+obsl)) == 0;
            obs(:,1:x-1) = false;
            obs(:,x+L*cosd(aoa)+1:l) = false;
            

            % figure;
            % imagesc(obs); colorbar;
            % axis equal;
        end
    end

    methods (Access=private)
        function calculateFunctions(obj, aoa, L)
            arguments
                obj
                aoa double; % angle of attack (degrees)
                L double;
            end

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

            xu = @(t) L*(t - yt(t).*dyc(t)./sqrt(1+dyc(t).^2))*cosd(aoa) + L*(yc(t) + yt(t)./(1+dyc(t).^2).^.5)*sind(aoa);
            dxu = @(t) L*(1 - (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*cosd(aoa) ...
                       + L*(dyc(t) + (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*sind(aoa);
            yu = @(t) L*(yc(t) + yt(t)./(1+dyc(t).^2).^.5)*cosd(aoa) - L*(t - yt(t).*dyc(t)./sqrt(1+dyc(t).^2))*sind(aoa);
            dyu = @(t) L*(dyc(t) + (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*cosd(aoa) ...
                       - L*(1 - (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*sind(aoa);
            xl = @(t) L*(t + yt(t).*dyc(t)/sqrt(1+dyc(t).^2))*cosd(aoa) + L*(yc(t) - yt(t)./(1+dyc(t).^2).^.5)*sind(aoa);
            dxl = @(t) L*(1 + (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*cosd(aoa) ...
                       + L*(dyc(t) - (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*sind(aoa);
            yl = @(t) L*(yc(t) - yt(t)./(1+dyc(t).^2).^.5)*cosd(aoa) - L*(t + yt(t).*dyc(t)/sqrt(1+dyc(t).^2))*sind(aoa);
            dyl = @(t) L*(dyc(t) - (dyt(t)./(1+dyc(t).^2).^.5 - yt(t).*(dyc(t).*d2yc(t))/(1+dyc(t).^2).^1.5))*cosd(aoa) ...
                       - L*(1 + (dyt(t).*dyc(t)./(1+dyc(t).^2).^.5 + d2yc(t)./(1+dyc(t).^2).^1.5))*sind(aoa);
            
            obj.wu = @(x,y) integral(@(t) ((xu(t)-x).*dyu(t) - (yu(t)-y).*dxu(t))./((xu(t)-x).^2 + (yu(t)-y).^2), 0,1, 'ArrayValued',true);
            obj.wl = @(x,y) integral(@(t) ((xl(t)-x).*dyl(t) - (yl(t)-y).*dxl(t))./((xl(t)-x).^2 + (yl(t)-y).^2), 0,1, 'ArrayValued',true);

        end
    %     function mustBeValid(name)
    %         if ~(name>=1000 && name <= 9999)
    %             error("invalid airfoil name. must be 4 digits")
    %         end
    %     end
    end
end