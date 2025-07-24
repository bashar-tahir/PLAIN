% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
%
%   Cite as:
%   B. Tahir, P. Svoboda, and M. Rupp, "PLAIN: Scalable Estimation
%   Architecture for Integrated Sensing and Communication," arXiv:
%   2503.21242, 2025. [Online]. Available: https://arxiv.org/abs/2503.21242
% -----------------------------------------------------------------------
%   The code implements the PLAIN architecture and demonstrates its use 
%   in an example scenario, where the sensing is carried out over the 
%   dimensions of space (azimuth), frequency, and time (OFDM symbols)
%
%%
clc; clear;

%% Parameters Setup
% Transmission
centerFreq    = 26e+09;                     % Center frequency
N_Antennas    = 16;                         % Number of receive antennas
N_Subcarriers = 180;                        % Number of subcarriers
N_TimeSymbols = 560;                        % Number of OFDM symbols
df            = 60e3;                       % Subcarrier spacing
FFT_Size      = N_Subcarriers*1.25;         % FFT size. Here we apply a bigger FFT compared to the number of subcarriers to virtually create fractional delays
samplingRate  = FFT_Size*df;                % Resulting sampling rate
T_Sample      = 1/samplingRate;             % Resulting sampling time
CPlength      = 1/(df*14);                  % Cyclic-prefix length
dt            = (1/df) + CPlength;          % Resulting time symbol spacing
c             = physconst('LightSpeed');    % Speed of light
wavelength    = c/centerFreq;               % Resulting wavelength
da            = wavelength/2;               % Chosen antenna spacing

% Environment 
N_P           = 6;                          % Number of objects/paths
angularRange  = [30, 150];                  % Angular range (deg) over which the N_P objects are uniformly spaced
distanceRange = [50, 400];                  % Distance range (m) over which the N_P objects are uniformly spaced
velocityRange = [0, 25];                    % Velocity (km/h) range over which the N_P objects are uniformly spaced (this is modified the Channel.m object)

% Simulation sweep
SNRvec_dB = -20:5:20;                       % SNR range
nReps     = 20;                             % Number of simulation repetitions

%% PLAIN
compressionMethod = "Averaging";            % "None", "Decimation", "Averaging", "Decimation2", "Smoothing"
fusionMethod      = "TensorOMP";            % "TensorLS", "TensorOMP"
useTrueN_P        = 1;                      % 0 to never use true N_P, 1 to use true N_P only in the final selection, 2 to use true N_P in all stages (estimation + fusion & selection)
Deltas            = [1, 4, 14];             % Compression factors across [space, frequency, time]

t1 = tic;
%% Simulation
% OFDM and channel objects
Modulator  = OFDM(N_Subcarriers, N_TimeSymbols, df, FFT_Size, CPlength);
ChannelObj = Channel(N_Subcarriers, N_TimeSymbols, df, dt, centerFreq, samplingRate, Modulator.N_SamplesTotal, N_P, [1, N_Antennas], angularRange, distanceRange, velocityRange);

% Generate a test realization and print out the resulting parameters
ChannelObj.generateChannel();
disp(['Scenario setup (N_P = ' num2str(N_P) '):'])
disp(num2str([ChannelObj.Paths.AOAd.'; ChannelObj.Paths.Distances; ChannelObj.Paths.Vels]))

% Sweep and results parameters
Pn = physconst('Boltzmann') * 296 * samplingRate * (N_Subcarriers/FFT_Size);    % Noise power
PtVec = 10.^(SNRvec_dB/10) * Pn / sum(10.^(-ChannelObj.Paths.PathLosses/10));   % Transmit power vector
lenPtVec = length(PtVec);
MSE = zeros(length(PtVec), nReps, N_P, 3);
MSE_CRB = zeros(length(PtVec), nReps, N_P, 3);

%% Repetition loop (use parfor here)
for iRep = 1:nReps
    rng(iRep, 'twister');
    % Generate a new channel realization
    t2 = tic;
    ChannelObj.generateChannel();
    r2 = toc(t2);
    % To save time, we apply the same channel realization to different SNR points
    for iPt = 1:lenPtVec
        rng(iRep, 'twister');
        Pt = PtVec(iPt);  

        % Generate transmit signal
        dataTimeFreq = reshape(randi([0,1], N_Subcarriers*N_TimeSymbols, 1)*2-1, N_Subcarriers, N_TimeSymbols); % Random BPSK symbols
        UEtimeSignal = sqrt(Pt) * Modulator.Modulate(dataTimeFreq);
        
        % Propagation
        yTotal = cell(N_Antennas, 1);
        for r = 1:N_Antennas
            % Convolve with the channel
            y = transpose(ChannelObj.ConvMatrix{r}*UEtimeSignal);
            
            % Add noise
            y = y + sqrt(Pn/2)*complex(randn(size(y)), randn(size(y)));
            
            % Demodulation
            y = y(ChannelObj.InitialDelayIdx:end);
            yTotal{r} = Modulator.Demodulate(y);
        end
        
        % Equalization and channel estimation
        H_total = zeros(1,N_Antennas,N_Subcarriers,N_TimeSymbols); % First dimension corresponds to the snapshots. Here, it is 1.
        for r = 1:N_Antennas
            H_est = yTotal{r} ./ dataTimeFreq;
            
            % Correct for time offset
            kVec = (0:N_Subcarriers-1).';
            H_est = H_est .* exp(-1i*2*pi*ChannelObj.InitialDelay*df*kVec);

            H_total(1,r,:,:) = H_est;
        end
        % The resulting H_total is of dim (S x N_1 x N_2 x ... x N_M), where
        % S is the number of snapshots. In this example, we have S = 1
        
        %% Parameters Estimation
        % PLAIN
        PLAINObj = PLAIN(compressionMethod,fusionMethod, useTrueN_P, N_P, centerFreq, da, df, dt, Deltas);
        t3 = tic;
        [N_P_est, EstiAOA, EstiRange, EstiVel, EstiGains, H_clean] = PLAINObj.extractObjects(H_total);
        r3 = toc(t3);
        
        %% RMSE and CRB Calculation       
        % Sort based on angle
        [~, I1] = sort(EstiAOA);
        [~, I2] = sort(ChannelObj.Paths.AOAd);
        % [~, I1] = sort(EstiRange);
        % [~, I2] = sort(Distance);
        
        % Store the MSE results, if length matches
        if N_P_est == N_P
            MSE_SimAOA = abs(EstiAOA(I1) - ChannelObj.Paths.AOAd(I2)).^2;
            MSE_SimDis = abs(EstiRange(I1)*1000 - ChannelObj.Paths.Distances(I2).').^2;
            MSE_SimVel = abs(EstiVel(I1)*1000/3600 - ChannelObj.Paths.Vels(I2).'*1000/3600).^2; 
            MSE(iPt, iRep, :) = [MSE_SimAOA; MSE_SimDis; MSE_SimVel];
        else
            MSE(iPt, iRep, :) = [nan(N_P, 1); nan(N_P, 1); nan(N_P, 1)];
        end
        
        t4 = tic;  
        % CRB
        % Construct the response matrices and their derivatives
        A = exp(1i*2*pi*cos(ChannelObj.Paths.AOAd.'*pi/180)*da.*(0:(N_Antennas-1)).'/wavelength);
        D = exp(-1i*2*pi*(ChannelObj.Paths.Distances*1000/c)*df.*(0:(N_Subcarriers-1)).');
        V = exp(1i*2*pi*(centerFreq*ChannelObj.Paths.Vels*1000/3600/c)*dt.*(0:(N_TimeSymbols-1)).');
        DA = (-1i*2*pi*pi/180*sin(ChannelObj.Paths.AOAd.'*pi/180)*da.*(0:(N_Antennas-1)).'/wavelength) .* A;
        DF = (-1i*2*pi*(ones(size(ChannelObj.Paths.Distances))/c)*df.*(0:(N_Subcarriers-1)).') .* D;
        DT = (1i*2*pi*(centerFreq*ones(size(ChannelObj.Paths.Vels))/c).*dt.*(0:(N_TimeSymbols-1)).') .* V;
        % Use the tensor structure to simplify the calculation
        CRB_AOA = diag(calcCRBTensor(1, Pt*ChannelObj.P, {DA, D, V, A}, Pn));
        CRB_Dis = diag(calcCRBTensor(1, Pt*ChannelObj.P, {DF, A, V, D}, Pn));
        CRB_Vel = diag(calcCRBTensor(1, Pt*ChannelObj.P, {DT, A, D, V}, Pn));
        MSE_CRB(iPt, iRep, :) = [CRB_AOA; CRB_Dis; CRB_Vel];
        r4 = toc(t4);
        
        % Simulation display for current repetition and sweep
        disp(['---------' ' Rep: ' num2str(iRep) ', SNR = ' num2str(10 * log10(Pt .* sum(10.^(-ChannelObj.Paths.PathLosses/10)) / Pn))  ' dB ---------'])
        disp(['Channel generation: ' num2str(r2) ' seconds.']);
        disp(['Algorithm: ' num2str(r3) ' seconds.']);
        disp(['CRB: ' num2str(r4) ' seconds.']);
        disp(' ')
    end
end
disp(['Total time: ' num2str(toc(t1)) ' seconds.']);

%% Plots
figure(1)
tiled = tiledlayout('horizontal');
% tiled.TileSpacing = 'compact';
tiled.Padding = 'compact';

% AOA RMSE
nexttile;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE(:,:,:,1), 3), 2, "omitnan")).','-ob', "LineWidth", 1); hold on;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE_CRB(:,:,:,1), 3), 2, "omitnan")).','--b', "LineWidth", 1); hold off;
ylabel('AOA RMSE [deg]')
xlabel('SNR [dB]')
grid on;
legend(["PLAIN", "CRB"]);

% Distance RMSE
nexttile;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE(:,:,:,2), 3), 2, "omitnan")).','-ob', "LineWidth", 1); hold on;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE_CRB(:,:,:,2), 3), 2, "omitnan")).','--b', "LineWidth", 1); hold off;
ylabel('Distance RMSE [m]')
xlabel('SNR [dB]')
grid on;

% Velocity RMSE
nexttile;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE(:,:,:,3), 3), 2, "omitnan")).','-ob', "LineWidth", 1); hold on;
semilogy(SNRvec_dB, sqrt(mean(mean(MSE_CRB(:,:,:,3), 3), 2, "omitnan")).','--b', "LineWidth", 1); hold off;
ylabel('Velocity RMSE [m/s]')
xlabel('SNR [dB]')
grid on;
