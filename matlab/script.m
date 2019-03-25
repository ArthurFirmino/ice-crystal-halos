%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
n_rays = 16e9;
n_threads = 12;
crystal = Crystal.Plate();
% crystal = Crystal.Column();

fprintf('Sampling %d rays...\n', n_rays);
tic
pf_parameters = [256,256,256];
phase_function = ics_function_mex(n_rays, n_threads, crystal, pi/8, true, pf_parameters);
% phase_function = ics_function_mex(n_rays, n_threads, crystal, pi, false, pf_parameters);
writeout_mex(phase_function, 'data/tabulated_pf.data');
toc