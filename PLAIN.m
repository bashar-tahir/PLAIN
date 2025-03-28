%% ------------------------------------------------------------------------
%   (c) 2025 Bashar Tahir, bashar.tahir@tuwien.ac.at
%   Institute of Telecommunications, TU Wien
%   https://www.tuwien.at/etit/tc/en/
% -----------------------------------------------------------------------
%   Implements the PLAIN architecture
%
%%
classdef PLAIN
    % Comp_ressed Decoupl_ed Estima_tion and I_nput-Based Fusion_ (PLAIN)
    
    properties
        CenterFrequency         % Center frequency
        da                      % Antenna spacing
        df                      % Subcarrier spacing
        dt                      % Time symbol spacing
        useTrueN_P              % 0 to never use true NP, 1 to use true NP only in the final selection, 2 to use true NP in all stages (estimation + fusion & selection)
        N_P                     % True number of objects
        CompressionMethod       % Compression method ("None", "Decimation", "Averaging", "Decimation2", "Smoothing")
        FusionMethod            % Fusion method ("TensorLS", "TensorOMP")
        Deltas                  % Compression factors
    end
    
    methods
        function obj = PLAIN(compressionMethod, fusionMethod, useTrueN_Pin, N_Pin, centerFrequency, dain, dfin, dtin, Deltasin)
            obj.CompressionMethod   = compressionMethod;
            obj.FusionMethod        = fusionMethod;
            obj.useTrueN_P          = useTrueN_Pin;
            obj.N_P                 = N_Pin;
            obj.CenterFrequency     = centerFrequency;
            obj.da                  = dain;
            obj.df                  = dfin;
            obj.dt                  = dtin;
            obj.Deltas              = Deltasin;
        end
        
        function [N_P_est, EstiAOA, EstiRange, EstiVel, EstiGains, H_clean] = extractObjects(obj, H_total)
            % Input H_total is of dim (S x N_1 x N_2 x ... x N_M), where
            % S is the number of snapshots. Snapshots are the first dim
            % inorder to avoid array collapse in MATLAB when the last dim
            % is singular (in the case of single snapshot).
            %
            % Output is the set of paired estimated angles (deg), distances (km),
            % Velocities (km/h) and path gains, as well as a cleaned-up H
            %
            c          = physconst('LightSpeed');
            wavelength = c/obj.CenterFrequency;

            %% Compression
            H = obj.compress(H_total, obj.CompressionMethod, [1, obj.Deltas]);

            %% Estimation
            [N_S, cN_R, cN_Subcarriers, cN_TimeSymbols] = size(H);

            if obj.CompressionMethod == "Smoothing"
                obj.Deltas(1) = 1;
                obj.Deltas(2) = 1;
                obj.Deltas(3) = 1;
            end
            
            % Estimation via root-MUSIC and AIC
            % Angle
            Hm = tunfold(H, 2);
            [cand_Angles, ~, ~] = ODE_MUSIC(Hm, obj.useTrueN_P, obj.N_P);
            cand_Angles = acos(cand_Angles*wavelength/(2*pi*obj.da*obj.Deltas(1)))*180/pi;      % radial range from -pi to pi
            
            % Distance
            Hm = tunfold(H, 3);
            [tempAngles, ~, ~] = ODE_MUSIC(Hm, obj.useTrueN_P, obj.N_P);
            tempAngles(tempAngles>0) = -pi - (pi - tempAngles(tempAngles>0));                   % radial range from 0 to 2*pi
            cand_Ranges = -tempAngles*c/(2*pi*obj.df*obj.Deltas(2)*1000);                  
            
            % Velocity
            Hm = tunfold(H, 4);
            [cand_Vels, ~, ~] = ODE_MUSIC(Hm, obj.useTrueN_P, obj.N_P);
            cand_Vels = cand_Vels / (2*pi*obj.Deltas(3)*obj.dt*obj.CenterFrequency*1000/3600/c); % radial range from -pi to pi

            % Reconstruct response matrices
            A = exp(1i*2*pi*cos(cand_Angles.'*pi/180)*obj.Deltas(1)*obj.da.*(0:(cN_R-1)).'/wavelength);
            D = exp(-1i*2*pi*(cand_Ranges.'*1000/c)*obj.Deltas(2)*obj.df.*(0:(cN_Subcarriers-1)).');
            V = exp(1i*2*pi*(obj.CenterFrequency*cand_Vels.'*1000/3600/c)*obj.Deltas(3)*obj.dt.*(0:(cN_TimeSymbols-1)).');
              
            %% Fusion and Selection
            if obj.useTrueN_P == 0
                N_P_est = max([length(cand_Angles), length(cand_Ranges), length(cand_Vels)]);
            else
                N_P_est = obj.N_P;
            end
            
            switch obj.FusionMethod
                case 'TensorLS'
                    B     = tmul(H, {eye(size(H,1)), pinv(A), pinv(D), pinv(V)});
                    meanB = mean(abs(B).^2, 1);
                    [~, idx] = sort(meanB(:), 'descend');
                    N_P_est = min(length(idx), N_P_est);
                    [~, a, f, t] = ind2sub(size(meanB), idx(1:N_P_est));
                    EstiGains = sqrt(meanB(idx(1:N_P_est)));
                    % Obtain a cleaned up version, by removing the so far detected objects
                    H_clean = [];
                    % H_clean = H - tmul(B, {eye(size(H,1)), A(:, unique(a)), D(:, unique(f)), V(:, unique(t))});
                case 'TensorOMP'
                    [a, f, t, EstiGains, H_clean] = TensorOMP(H, {A, D, V}, N_P_est);
                otherwise
                    error('Pairing method not supported.');
            end
           
            EstiRange = cand_Ranges(f);
            EstiVel   = cand_Vels(t);
            EstiAOA   = cand_Angles(a);
        end
    end
    methods(Static)
        %%
        function H_c = compress(H, method, Deltas)
            % Compress the input H according a certain method and the
            % compression factors in Deltas
            [N_S, N_R, N_Subcarriers, N_TimeSymbols] = size(H);
            sizeArray = size(H);
            switch method
                case 'None'
                    H_c = H;
                case 'Decimation'
                    H_c = H(:, 1:Deltas(2):end, 1:Deltas(3):end, 1:Deltas(4):end);
                case 'Averaging'
                    Jcell = cell(1,length(sizeArray));
                    for m = 1:length(sizeArray)
                        J = zeros(sizeArray(m), sizeArray(m)/Deltas(m));
                        for i = 0:(sizeArray(m)/Deltas(m)-1)
                            e = zeros(sizeArray(m), 1);
                            e(i*Deltas(m)+1:((i+1)*Deltas(m))) = 1;
                            J(:,i+1) = e;
                        end
                        Jcell{m} = J.' ./ Deltas(m);
                    end
                    H_c = tmul(H, Jcell);
                case 'Decimation2'
                    % Limit number of shifts to 100
                    N_VS = min(prod(Deltas), 100);
                    H_c = zeros([N_VS sizeArray(2:end)./Deltas(2:end)]);
                    [~, I2, I3, I4] = ind2sub(Deltas, randperm(prod(Deltas), N_VS));
                    for s = 1:N_VS
                        H_c(s,:,:,:) = H(:, I2(s):Deltas(2):end, I3(s):Deltas(3):end, I4(s):Deltas(4):end);
                    end
                case 'Smoothing'
                    N_RS           = N_R / Deltas(2);
                    N_SubcarriersS = N_Subcarriers / Deltas(3);
                    N_TimeSymbolsS = N_TimeSymbols / Deltas(4);

                    % Skip factors for uniform sampling of the shifts. For
                    % the parameter settings in the paper, the following
                    % jumps are applied to produce a total of 100 snapshots
                    sA = 1;
                    sF = 10;
                    sT = sF*3;

                    N_S = length(1:sA:(N_R-N_RS+1))*length(1:sF:(N_Subcarriers-N_SubcarriersS+1))*length(1:sT:(N_TimeSymbols-N_TimeSymbolsS+1));
                    H_c = zeros(N_S, N_RS, N_SubcarriersS, N_TimeSymbolsS);
                    k = 1;
                    for iR = 1:sA:(N_R-N_RS+1)
                        for iF = 1:sF:(N_Subcarriers-N_SubcarriersS+1)
                            for iT = 1:sT:(N_TimeSymbols-N_TimeSymbolsS+1)
                                H_c(k,:,:,:) = H(1,iR:(iR+N_RS-1), iF:(iF+N_SubcarriersS-1), iT:(iT+N_TimeSymbolsS-1));
                                k = k + 1;
                            end
                        end
                    end
                otherwise
                    error('Compression type not known.')
            end
        end
    end
end

