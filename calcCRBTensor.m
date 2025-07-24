%% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
% -----------------------------------------------------------------------
%   Approximate deterministic CRB calculation using the tensor structure
%   for three factor matrices.
%   Following Equation (8.110) of "H. Van Trees, Optimum Array Processing: 
%   Part IV of Detection, Estimation, and Modulation Theory. Wiley, 2004."
%
%%
function CRB = calcCRBTensor(N_s, Sf, Z, sigma2)
    D  = Z{1};
    A1 = Z{2};
    A2 = Z{3};
    A3 = Z{4};
    DhD = (D'*D) .* (A1'*A1) .* (A2'*A2);
    DhV = (D'*A3) .* (A1'*A1) .* (A2'*A2);
    H = DhD - DhV * ((A1'*A1) .* (A2'*A2) .* (A3'*A3))^-1 * DhV';
    CRB = (real(Sf .* H.')^-1) .* sigma2 / (2 * N_s);
end