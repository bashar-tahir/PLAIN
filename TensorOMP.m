%% ------------------------------------------------------------------------
%   (c) 2025 Bashar Tahir, bashar.tahir@tuwien.ac.at
%   Institute of Telecommunications, TU Wien
%   https://www.tuwien.at/etit/tc/en/
% -----------------------------------------------------------------------
%   Tensor-OMP fusion for four dimensions (incl. snapshots dimension)
%
%%
function [I1, I2, I3, betas, H_clean] = TensorOMP(H, U, N_P_est)
    R = H;
    I1 = [];
    I2 = [];
    I3 = [];
    
    N_S = size(H,1); % Number of snapshots
    A = U{1};
    D = U{2};
    V = U{3};
    for iter = 1:N_P_est
        % Find strongest path
        X  = tmul(R, {eye(N_S), A', D', V'});
        meanX = sum(abs(X).^2, 1);
        [~, idx] = max(meanX(:)); 
        [~, iA, iF, iT] = ind2sub(size(meanX), idx);
        I1 = [I1; iA];
        I2 = [I2; iF];
        I3 = [I3; iT];
        I1u = unique(I1);
        I2u = unique(I2);
        I3u = unique(I3);

        % Estimate the core coefficients
        B = tmul(H, {eye(N_S), pinv(A(:, I1u)), pinv(D(:, I2u)), pinv(V(:, I3u))});

        % Construct an approxmiation of H
        H_approx = tmul(B, {eye(N_S), A(:, I1u), D(:, I2u), V(:, I3u)});

        % Update the residual
        R = H - H_approx;
    end
    H_clean = R;

    % Path gains
    meanB = mean(abs(B).^2, 1);
    [~, idx] = sort(meanB(:), 'descend');
    betas = sqrt(meanB(idx(1:N_P_est)));
end

