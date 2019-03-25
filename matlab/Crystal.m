%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
classdef Crystal
    % CRYSTAL class, represents an ice crystal
    % Contains crystal geometry and performs ray sampling and tracing
    
    properties
        BoundingSphereR
        n_r
        TriangleList
    end
    
    properties(Constant)
        % wavelengths(um) and corresponding indices of refraction
        wl_array = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, ...
            0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90];
        n_r_array = [1.3339, 1.3249, 1.3194, 1.3157, 1.3130, 1.3110, ...
            1.3094, 1.3080, 1.3069, 1.3059, 1.3049, 1.3040, 1.3032];
        wavelength_min = 0.36; % Keep these values consistent with
        wavelength_max = 0.83; % TabulatedPF
    end
    
    methods(Static)
        function obj = HexagonalCylinder(L, W)
            % Create Crystal object, set properties
            obj = Crystal;
            obj.BoundingSphereR = sqrt((L/2)^2 + (W/2)^2);
            obj.n_r = 1.309;
            obj.TriangleList = zeros(3,3,24);
            % Define triangles
            for i = 0:1:5
                % Angles
                b = i*2*pi/6;
                c = (i+1)*2*pi/6;
                % Top and bottom triangles
                A = [0,0,L/2];
                B = [cos(b),sin(b),L/2].*[W/2,W/2,1];
                C = [cos(c),sin(c),L/2].*[W/2,W/2,1];
                obj.TriangleList(:,:,1+4*i) = [A;B;C];
                A = A.*[1,1,-1];
                B = B.*[1,1,-1];
                C = C.*[1,1,-1];
                obj.TriangleList(:,:,2+4*i) = [A;C;B];
                % 2 side triangles
                A = [cos(b),sin(b),L/2].*[W/2,W/2,-1];
                B = [cos(c),sin(c),L/2].*[W/2,W/2,1];
                C = [cos(b),sin(b),L/2].*[W/2,W/2,1];
                obj.TriangleList(:,:,3+4*i) = [A;B;C];
                C = [cos(c),sin(c),L/2].*[W/2,W/2,-1];
                obj.TriangleList(:,:,4+4*i) = [A;C;B];
            end
        end
        function obj = Plate()
            % Plate dimensions
            L = 35.9; W = 175.1;
            obj = Crystal.HexagonalCylinder(L, W);
        end
        function obj = Column()
            % Plate dimensions
            L = 240.6; W = 108.2;
            obj = Crystal.HexagonalCylinder(L, W);
        end
    end
    
    methods
        
        function [t_min,i_min] = closestHit(obj, R)
            n = size(obj.TriangleList,3);
            t = Inf*ones(n,1);
            rO = R.o;
            rV = R.d;
            eps = 0.0000001;
            % Perform ray triangle intersection
            for i = 1:n
                A = obj.TriangleList(1,:,i);
                B = obj.TriangleList(2,:,i);
                C = obj.TriangleList(3,:,i);
                % Moller-Trumbore algorithm
                edge1 = B-A;
                edge2 = C-A;
                h = cross(rV,edge2);
                a = edge1*(h');
                if a > -eps && a < eps
                    continue
                end
                f = 1/a;
                s = rO - A;
                u = f*(s*(h'));
                if u < 0.0 || u > 1.0
                    continue
                end
                q = cross(s,edge1);
                v = f*(rV*(q'));
                if v < 0.0 || u + v > 1.0
                    continue
                end
                T = f*(edge2*(q'));
                if T > eps
                   t(i) = T; 
                end
            end
            % Find closest hit, if any
            t_min = Inf;
            i_min = 0;
            for i = 1:n
                if t(i) < t_min
                    t_min = t(i);
                    i_min = i;
                end
            end
        end
        
        function [R,hit] = simulateRay(obj, R)
            obj.n_r = interp1(obj.wl_array, obj.n_r_array, R.wl);
            [t,i] = obj.closestHit(R);
            if( t<Inf )
                hit = true;
            else
                hit = false;
            end
            while(t<Inf)
                % Calculate normal
                N = cross( obj.TriangleList(2,:,i)-obj.TriangleList(1,:,i),...
                obj.TriangleList(3,:,i)-obj.TriangleList(1,:,i));
                N = (1-2*R.in_object)*N / norm(N);
                % Solve Snells Law, and Fresnel equations
                n_i = R.n_r;
                n_t = (1-R.in_object)*obj.n_r + R.in_object*1.0;
                theta_i = atan2(norm(cross(N,-R.d)),dot(N,-R.d));
                sintheta_t = (n_i/n_t)*sin(theta_i);
                theta_t = 0.0;
                T = -0.0000001;
                if sintheta_t <= 1.0
                    theta_t = asin((n_i/n_t)*sin(theta_i));
                    % Solve for T (Fresnel equations)
                    rpl = (n_t*cos(theta_i) - n_i*cos(theta_t))/...
                        (n_t*cos(theta_i) + n_i*cos(theta_t));
                    rpr = (n_i*cos(theta_i) - n_t*cos(theta_t))/...
                        (n_i*cos(theta_i) + n_t*cos(theta_t));
                    T = 1.0 - 0.5*(rpl*rpl + rpr*rpr);
                end
                % Russian roulette
                offset = 0.0000001;
                if rand < T
                    R.o = R.o + (t+offset)*R.d;
                    R.d = (n_i/n_t)*R.d + ((cos(theta_i)*n_i/n_t)-...
                        sqrt(1-sin(theta_t)*sin(theta_t)))*N;
                    R.n_r = n_t;
                    R.in_object = ~R.in_object;
                else
                    R.o = R.o + (t-offset)*R.d;
                    R = R.reflect(N);
                end
                [t,i] = obj.closestHit(R);
            end
        end
        
        function [w_l, theta_i, theta_s, delta_phi, hit] = sample_incident_ray(obj, M)
            % Sample incident ray angles
            theta_i = acos(2*rand -1);
            phi_i = 2*pi*rand;
            % Add small random rotations to crystal
            M_inv = M';
            % Calculate ray direction incident ray direction
            dir = [cos(phi_i)*sin(theta_i),sin(phi_i)*sin(theta_i),cos(theta_i)];
            dir = (M_inv*dir')';
            org = [0,0,0]+obj.BoundingSphereR*dir;
            % Offset origin by sampling project bounding sphere
            tangent1 = [cos(phi_i)*cos(theta_i),sin(phi_i)*cos(theta_i),-sin(theta_i)];
            tangent2 = [-sin(phi_i),cos(phi_i),0];
            tangent1 = (M_inv*tangent1')';
            tangent2 = (M_inv*tangent2')';

            radius_offset = obj.BoundingSphereR * sqrt(rand);
            angle_offset = 2*pi*rand;
            offset_vector = radius_offset*(cos(angle_offset)*tangent1 +...
                                        sin(angle_offset)*tangent2);
            org = org + offset_vector;
            
            % Sample wavelength from range [.38,.75]
            w_l = obj.wavelength_min + rand*(obj.wavelength_max-obj.wavelength_min);
            R_i = Ray(org,-dir,1.0,w_l,false);
            % Simulate ray until it exits bounding sphere
            if true %t < Inf
               [R_s,hit] = obj.simulateRay(R_i); 
            end
            R_s.d = (M*R_s.d')';
            % Record scattered direction
            [phi_s, theta_s, ~] = cart2sph(R_s.d(1),R_s.d(2),R_s.d(3));
            theta_s = -theta_s + pi/2;
            delta_phi = angdiff(phi_i,phi_s);
        end
        
        function visualize(obj)
            % Plot crystal
            hold on
            for i = 1:1:24
                X = obj.TriangleList(:,1,i);
                Y = obj.TriangleList(:,2,i);
                Z = obj.TriangleList(:,3,i);
                plot3([X;X(1)],[Y;Y(1)],[Z;Z(1)],'k-')
                hold on
            end
            axis([-130 130 -130 130 -130 130])
            xlabel('x (\mum)')
            ylabel('y (\mum)')
            zlabel('z (\mum)')
            hold off
        end
    end
end