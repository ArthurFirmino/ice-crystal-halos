%{
    This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

    Copyright (c) 2019 by Arthur Firmino.
%}
classdef TabulatedPF
    % TABULATEDPF Used to integrate simulated light rays into a phase
    % function, and write to file
    properties
        nThetaI, nThetaO, nDeltaPhi
        phaseFunctionArray
        nRaysArray
        nRaysHitArray
        spectralPhaseFunctionArray
    end
    
    methods(Static)
        function obj = TabulatedPF(nThetaI, nThetaO, nDeltaPhi)
            C = Constants;

            obj.nThetaI = nThetaI;
            obj.nThetaO = nThetaO;
            obj.nDeltaPhi = nDeltaPhi;
            obj.phaseFunctionArray = zeros(nThetaI, nThetaO, nDeltaPhi);
            obj.nRaysArray = zeros(nThetaI, 1);
            obj.nRaysHitArray = zeros(nThetaI, 1);
            if( C.SPECTRUM_SAMPLES > 1)
                obj.spectralPhaseFunctionArray = zeros(nThetaI, nThetaO, nDeltaPhi, ...
                    C.SPECTRUM_SAMPLES);
            end
        end
    end
    
    methods        
        function obj = addsample(obj, sampleInfo)
            C = Constants;

            % Convert sample_info data to ranges [0,1]
            ti = 0.5 * (1 + cos( sampleInfo.thetaI ));
            to = 0.5 * (1 + cos( sampleInfo.thetaO ));
            dp = abs( sampleInfo.deltaPhi ) / pi;
            wl = ( sampleInfo.wavelength - C.WAVELENGTH_MIN ) / ...
                ( C.WAVELENGTH_MAX - C.WAVELENGTH_MIN );
            
            % Compute left indices, 0-based
            iIndex = min(floor(ti * ( obj.nThetaI - 1 )), obj.nThetaI - 2);
            jIndex = min(floor(to * ( obj.nThetaO - 1 )), obj.nThetaO - 2);
            kIndex = min(floor(dp * ( obj.nDeltaPhi - 1 )), obj.nDeltaPhi - 2);
            lIndex = 1;
            if( C.SPECTRUM_SAMPLES > 1 )
                lIndex = min(floor(wl * C.SPECTRUM_SAMPLES), ...
                    C.SPECTRUM_SAMPLES - 1) + 1;
                lIndex = 1 + C.SPECTRUM_SAMPLES - lIndex;
                if( C.SPECTRUM_SAMPLES == 3 && lIndex == 1)
                    lIndex = 3;
                else if( C.SPECTRUM_SAMPLES == 3 && lIndex == 3)
                    lIndex = 1;
                end
            end

            % Compute linear grid point contributions
            % Add them to corresponding arrays
            w0 = 0; w1 = 0; w2 = 0;
            for i = 1:2
                if( i == 1 ) w0 = 1 + double(iIndex) - (double(obj.nThetaI) - 1)*ti;
                else w0 = 1 - w0; end
                obj.nRaysArray(iIndex+i) = obj.nRaysArray(iIndex+i) + w0;
                if( sampleInfo.scattered )
                    obj.nRaysHitArray(iIndex+i) = obj.nRaysHitArray(iIndex+i) + w0;
                else
                    continue
                end
                for j = 1:2
                    if( j == 1 ) w1 = w0*(1 + double(jIndex) - (double(obj.nThetaO) - 1)*to);
                    else w1 = w0 - w1; end
                    for k = 1:2
                        if( k == 1 ) w2 = w1*(1 + double(kIndex) - (double(obj.nDeltaPhi) - 1)*dp);
                        else w2 = w1 - w2; end
                        edgeDouble = 1;
                        if( iIndex+i==1 || jIndex+j==1 || kIndex+k == 1 || ...
                            iIndex+i==obj.nThetaI || jIndex+j==obj.nThetaO || kIndex+k==obj.nDeltaPhi)
                            edgeDouble = 2;
                        end
                        obj.phaseFunctionArray(iIndex+i, jIndex+j, kIndex+k) = ...
                            obj.phaseFunctionArray(iIndex+i, jIndex+j, kIndex+k) + ...
                            edgeDouble * w2;
                        if( C.SPECTRUM_SAMPLES > 1) 
                            obj.spectralPhaseFunctionArray(iIndex+i, jIndex+j, kIndex+k, lIndex) = ...
                                obj.spectralPhaseFunctionArray(iIndex+i, jIndex+j, kIndex+k, lIndex) + ...
                                edgeDouble * w2;
                        end
                    end
                end
            end      
        end
        
        function obj = addphasefunction(obj, pf)
            C = Constants;
            
            obj.phaseFunctionArray = obj.phaseFunctionArray + ...
                pf.phaseFunctionArray;
            obj.nRaysArray = obj.nRaysArray + ...
                pf.nRaysArray;
            obj.nRaysHitArray = obj.nRaysHitArray + ...
                pf.nRaysHitArray;
            if(C.SPECTRUM_SAMPLES > 1)
                obj.spectralPhaseFunctionArray = ...
                    obj.spectralPhaseFunctionArray + ...
                    pf.spectralPhaseFunctionArray;
            end
        end
        
        function writetofile(obj, filename)
            C = Constants;
            
            for i=1:obj.nThetaI
               for j=1:obj.nThetaO
                   for k=1:obj.nDeltaPhi
                        obj.spectralPhaseFunctionArray(i,j,k,:) = ...
                        flip(obj.spectralPhaseFunctionArray(i,j,k,:));
                   end
               end
            end

            % Normalize pdf
            fprintf('Normalizing phase function...\n');
            dThetaO = 2 / (obj.nThetaO - 1);
            dDeltaPhi = 2 * pi / (obj.nDeltaPhi - 1);
            for l = 1:C.SPECTRUM_SAMPLES
                for i = 1:obj.nThetaI
                    intSpectral = 0;
                    intMono = 0;
                    for j = 1:obj.nThetaO
                        thetaOEdge = 1;
                        if( j == 1 || j == obj.nThetaO ) thetaOEdge = 0.5; end
                        for k = 1:obj.nDeltaPhi
                            deltaPhiEdge = 1;
                            if( k == 1 || k == obj.nDeltaPhi ) deltaPhiEdge = 0.5; end
                            intMono = intMono + dThetaO * thetaOEdge * dDeltaPhi * deltaPhiEdge * ...
                                obj.phaseFunctionArray(i,j,k);
                            if( C.SPECTRUM_SAMPLES > 1)
                                intSpectral = intSpectral + dThetaO * thetaOEdge * dDeltaPhi * deltaPhiEdge * ...
                                    obj.spectralPhaseFunctionArray(i, j, k, l);
                            end

                        end
                    end
                    if( intMono > 0 )
                        for j = 1:obj.nThetaO
                            for k = 1:obj.nDeltaPhi
                                obj.phaseFunctionArray(i,j,k) = ...
                                    obj.phaseFunctionArray(i,j,k) / intMono; 
                                if( C.SPECTRUM_SAMPLES > 1 && intSpectral > 0)
                                    obj.spectralPhaseFunctionArray(i, j, k, l) = ...
                                        obj.spectralPhaseFunctionArray(i, j, k, l) / intSpectral;
                                end
                            end
                        end
                    end
                end
            end

            % Compute directional sigma
            fprintf('Compute directional scattering coefficient...\n');
            sigmaScale = zeros(obj.nThetaI, 1);
            for i = 1:obj.nThetaI
                sigmaScale(i) = obj.nRaysHitArray(i) / obj.nRaysArray(i);
            end
            sigmaScale = sigmaScale / max(sigmaScale);

            % Write to file
            fprintf('Writing data to file...\n');
            
            fid = fopen(filename, 'w');
            if fid == -1
                error('Cannot open file for writing');
            end

            fwrite(fid, C.SPECTRUM_SAMPLES, 'uint32');
            fwrite(fid, obj.nThetaI, 'uint32');
            fwrite(fid, obj.nThetaO, 'uint32');
            fwrite(fid, obj.nDeltaPhi, 'uint32');
            
            pft = obj.phaseFunctionArray;
            pft = permute(pft, [3 2 1]);

            fwrite(fid, pft, 'float32');
            fwrite(fid, sigmaScale, 'float32');

            if( C.SPECTRUM_SAMPLES > 1) 
                pftS = obj.spectralPhaseFunctionArray;
                pftS = permute(pftS, [4 3 2 1]);
                fwrite(fid, pftS, 'float32');
            end
            fclose(fid);

            fprintf('Phase function succesfully written to file!\n');
        end
    end
end