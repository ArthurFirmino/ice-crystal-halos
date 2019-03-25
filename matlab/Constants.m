%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
classdef Constants
    properties(Constant)
        WAVELENGTH_MIN = 0.36;
        WAVELENGTH_MAX = 0.83;
        SPECTRUM_SAMPLES = 3;
    end
end

