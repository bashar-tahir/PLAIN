%% ------------------------------------------------------------------------
%   (c) 2025 Bashar Tahir, bashar.tahir@tuwien.ac.at
%   Institute of Telecommunications, TU Wien
%   https://www.tuwien.at/etit/tc/en/
% -----------------------------------------------------------------------
%   Sum-of-paths channel generation
%
%%
classdef Channel < handle
    % Channel object
    
    properties
        % Parameters
            N_Subcarriers       % Number of subcarriers
            N_TimeSymbols       % Number of OFDM symbols
            df                  % Subcarrier spacing
            dt                  % Time symbol spacing
            CenterFreq          % Center frequency
            Wavelength          % Wavelengt        
            SamplingRate        % Sampling rate
            N_Samples           % Total number of time samples
            N_P                 % Number of paths/objects
            AntennaConfig       % [N_TX, N_RX]
            P                   % Sources power autocorrelation matrix
            AngularRange        % Angular range over which the N_P objects are uniformly spaced
            DistanceRange       % Distance range over which the N_P objects are uniformly spaced
            VelocityRange       % Velocity range over which the N_P objects are uniformly spaced (this is modified below to have 50% static objects)
        % Output
            ConvMatrix      = []
            Paths           = []
            InitialDelayIdx = []
            InitialDelay    = []
    end
    
    methods
        function obj = Channel(varargin)
            obj.N_Subcarriers   = varargin{1};
            obj.N_TimeSymbols   = varargin{2};
            obj.df              = varargin{3};
            obj.dt              = varargin{4};
            obj.CenterFreq      = varargin{5};
            obj.Wavelength      = physconst('LightSpeed')/obj.CenterFreq;
            obj.SamplingRate    = varargin{6};
            obj.N_Samples       = varargin{7};
            obj.N_P             = varargin{8};
            obj.AntennaConfig   = varargin{9};
            obj.AngularRange    = varargin{10};
            obj.DistanceRange   = varargin{11};
            obj.VelocityRange   = varargin{12};
        end

        function PL = calcPathloss(obj,d)
            % 3GPP TR 38.901 version 15.0, UMa LOS, < dBPp, no SF corr
            PL = 28.0 + 22*log10(d) + 20*log10(obj.CenterFreq/1e9) + 4*0.74;
        end
        
        function generateChannel(obj)
            c = physconst('LightSpeed');
            T_Sample = 1/obj.SamplingRate;
            
            % Delays
            obj.Paths.Delays    = linspace(obj.DistanceRange(1), obj.DistanceRange(2), obj.N_P) / c;
            delayIdx            = round(obj.Paths.Delays/T_Sample)+1;
            obj.Paths.delayIdx  = delayIdx;
            obj.Paths.Delays    = T_Sample*(delayIdx-1); % Match to the possible taps of the channel (depends on the FFT size)

            % Gain coefficients
            obj.Paths.Distances   = obj.Paths.Delays * c;
            obj.Paths.PathLosses  = obj.calcPathloss(obj.Paths.Distances);
            phases                = exp(1i * rand(1, obj.N_P)*2*pi);
            RCS = 1;
            obj.Paths.betas      = sqrt(RCS).* phases .* sqrt(10.^(-obj.Paths.PathLosses/10));
                   
            % Doppler
            obj.Paths.Vels          = linspace(obj.VelocityRange(1), obj.VelocityRange(2), obj.N_P);
            obj.Paths.Vels(1:round(obj.N_P / (0.5*obj.N_P)):obj.N_P) = 0;                               % Set 50% of objects to static
            obj.Paths.DopplerFactor = (obj.Paths.Vels*1000/3600) / c;
            obj.Paths.Doppler       = obj.CenterFreq*obj.Paths.DopplerFactor;
            
            % Delay-Time
            minDelayIdx = min(delayIdx);
            obj.InitialDelayIdx = minDelayIdx;
            obj.InitialDelay    = (minDelayIdx-1)*T_Sample;

            % Space
            obj.Paths.AOAd  = linspace(obj.AngularRange(1), obj.AngularRange(2), obj.N_P).';
            obj.Paths.AOA   = pi * obj.Paths.AOAd / 180;
            da              = obj.Wavelength / 2;

            % Generate channel convolution matrix for each receive antenna r
            s1 = obj.N_Samples+minDelayIdx-1;
            s2 = obj.N_Samples;
            N_R = obj.AntennaConfig(2);
            obj.ConvMatrix = cell(N_R, 1);
            for r = 1:N_R
                obj.ConvMatrix{r} = sparse(s1, s2);
            end
            obj.P = RCS * diag(10.^(-obj.Paths.PathLosses/10));
            
            for p = 1:obj.N_P
                tracep = obj.Paths.betas(p) .* exp(1i*2*pi*obj.Paths.Doppler(p)*T_Sample*((delayIdx(p):s1)-1));
                for r = 1:N_R
                    obj.ConvMatrix{r} = obj.ConvMatrix{r} + sparse(delayIdx(p):s1, 1:obj.N_Samples-(delayIdx(p)-minDelayIdx), tracep .* exp(1i*2*pi*cos(obj.Paths.AOA(p))*da*(r-1)/obj.Wavelength), s1, s2);
                end
            end
        end

    end
end

