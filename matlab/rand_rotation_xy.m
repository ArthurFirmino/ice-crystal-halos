%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
function R = rand_rotation_xy(theta_range, normal)
    
    function M = rotation_matrix(u, theta)
        u = u / norm(u);
        ux = [0 -u(3) u(2) ; u(3) 0 -u(1) ; -u(2) u(1) 0 ];
        M = cos(theta)*eye(3) + sin(theta)*ux + (1-cos(theta))*kron(u,u');
    end

    phi_u_rand = 2*pi*rand;
    u_rand = [cos(phi_u_rand),sin(phi_u_rand),0];
    if normal
        theta_rand = theta_range*randn;
    else 
        theta_rand = theta_range*rand;
    end
    phi_z_rand = 2*pi*rand;
    R = rotation_matrix([0 0 1], phi_z_rand);
    R = rotation_matrix(u_rand, theta_rand) * R;
end

