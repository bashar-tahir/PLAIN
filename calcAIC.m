%% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
% -----------------------------------------------------------------------
%   Model-order selection based on AIC
%
%%
function N_P = calcAIC(lambda, N_samples)
    N_dim = length(lambda);
    IC = zeros(1, N_dim-1);
    for nA = 1:length(IC)
        nPara = 1 + nA * (2*N_dim - nA);
        sigman2_est = mean(lambda((nA+1):end));
        logf = N_samples * ((N_dim-nA)*log(sigman2_est) - sum(log(lambda((nA+1):end))));
        pen = nPara; 
        IC(nA) = logf + pen;
    end
    [~, N_P] = min(real(IC));
 
end

