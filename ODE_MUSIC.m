%% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
% -----------------------------------------------------------------------
%   One-dimensional estimation using root-MUSIC with FBA and AIC
%
%%
function [candidates, lambda, Vs] = ODE_MUSIC(Y, useTrueN_P, trueN_P)
    [N_dim, N_samples]  = size(Y);
    R = Y*Y'/N_samples;
    
    % Normalize for stability
    scaling = max(abs(R(:)));
    R = R ./ scaling;

    % Forward-backward averaging (FBA)
    J = fliplr(diag(ones(size(R,1),1)));
    R = 0.5*(R + J*conj(R)*J);

    [V, lambda] = eig(R, 'vector');
    [lambda, I] = sort(real(lambda), 'descend');
    V = V(:, I);
    
    % Determine N_P
    N_P = 1;
    if useTrueN_P < 2 
        N_P = calcAIC(lambda, N_samples);
    else
        N_P = trueN_P;
    end
    if N_P == 0; N_P = 1; end

    Vs = V(:, 1:N_P);
    
    % Root finding for off-grid parameter estimation
    zpoly = conv(Vs(:,1), conj(flipud(Vs(:,1))));
    for ii = 2:size(Vs,2)
        zpoly = zpoly + conv(Vs(:,ii), conj(flipud(Vs(:,ii))));
    end
    zpoly(N_dim) = zpoly(N_dim) - N_dim;
    rootz = roots(zpoly);
    rootz = rootz(abs(rootz)<=1);
    [~,I] = sort(abs(1-abs(rootz)));
    candidates = angle(rootz(I(1:N_P)));
end

