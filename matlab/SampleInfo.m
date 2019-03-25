%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
classdef SampleInfo
    % SAMPLEINFO, contains processed details of a ray sample
      
    properties
        thetaI, thetaO, deltaPhi
        wavelength
        scattered
    end
    
    methods
        function obj = SampleInfo(thetaI, thetaO, deltaPhi, ...
                wavelength, scattered)
            obj.thetaI = thetaI;
            obj.thetaO = thetaO;
            obj.deltaPhi = deltaPhi;
            obj.wavelength = wavelength;
            obj.scattered = scattered;
        end
    end
end

