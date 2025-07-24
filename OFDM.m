%% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
% -----------------------------------------------------------------------
%   Basic OFDM modulation and demodulation
%
%%
classdef OFDM < handle
    % OFDM Object
    
    properties
        % Parameters
            N_Subcarriers
            N_TimeSymbols
            df
            FFT_Size
            CP_Length
            N_CP_Samples
            SamplingRate
            SamplingTime
            N_SamplesTotal
    end
    
    methods
        function obj = OFDM(varargin)
            obj.N_Subcarriers   = varargin{1};
            obj.N_TimeSymbols   = varargin{2};
            obj.df              = varargin{3};
            obj.FFT_Size        = varargin{4};
            obj.CP_Length       = varargin{5};
            obj.SamplingRate    = obj.FFT_Size*obj.df;
            obj.SamplingTime    = 1/obj.SamplingRate;
            obj.N_CP_Samples    = round(obj.CP_Length / obj.SamplingTime);
            obj.N_SamplesTotal  = (obj.N_CP_Samples + obj.FFT_Size)*obj.N_TimeSymbols;
        end

        function timeSignal = Modulate(obj, inputTFGrid)
            baseGrid = zeros(obj.FFT_Size, obj.N_TimeSymbols);
            baseGrid(1:obj.N_Subcarriers, :) = inputTFGrid;
            baseGrid = ifft(baseGrid)*sqrt(obj.FFT_Size);
            outGrid  = [baseGrid(end-(obj.N_CP_Samples-1):end,:); baseGrid];
            timeSignal = outGrid(:);
        end

        function TFGrid = Demodulate(obj, inputTimeSignal)
            inGrid = reshape(inputTimeSignal, [], obj.N_TimeSymbols);
            TFGrid = fft(inGrid(obj.N_CP_Samples+1:end,:))*sqrt(1/obj.FFT_Size);
            TFGrid = TFGrid(1:obj.N_Subcarriers,:);
        end
    end
end

