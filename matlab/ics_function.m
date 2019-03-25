%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
function phase_function = ics_function(n_rays, n_threads, crystal, R_fun_input, R_fun_normal, pf_p)

    R_fun = @() rand_rotation_xy(R_fun_input, R_fun_normal);
    phase_function = TabulatedPF(pf_p(1), pf_p(2), pf_p(3));
    n_samples = 100;
    n_splits = floor(n_rays / (n_threads * n_samples));
    
    sampleArray = zeros(n_threads, n_samples,5);
    for i=1:n_splits
        parfor j=1:n_threads
            for k=1:n_samples
                M = R_fun();
                [w_l, theta_i, theta_s, delta_phi, hit] = crystal.sample_incident_ray(M);
                if( abs(theta_i-(pi-theta_s))<1e-6 && (abs(delta_phi)-pi)<1e-6 )
                    hit = false;
                end
                sampleArray(j,k,:) = [theta_i, theta_s, delta_phi, w_l, hit];
            end
        end
        for j=1:n_threads
            for k=1:n_samples
                s = SampleInfo(sampleArray(j,k,1),sampleArray(j,k,2),...
                    sampleArray(j,k,3),sampleArray(j,k,4),sampleArray(j,k,5));
                phase_function = phase_function.addsample(s);
            end
        end
    end
end

