%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
classdef Ray
    % RAY class, represents a ray used in CRYSTAL ray tracing
    
    properties
        o
        d
        n_r
        wl
        in_object
    end
    
    methods
        function obj = Ray(origin, direction, n_r, wl, in_object)
            obj.o = origin;
            obj.d = direction;
            obj.n_r = n_r;
            obj.wl = wl;
            obj.in_object = in_object;
        end
        
        function obj = reflect(obj, N)
           obj.d = obj.d - 2*dot(obj.d,N)*N;           
        end
        
        function visualize(obj, t)
           O = obj.o;
           D = obj.d;
           hold on
           plot3([O(1);O(1)+t*D(1)],[O(2);O(2)+t*D(2)],[O(3);O(3)+t*D(3)])
           hold off
        end
    end
end