%% ------------------------------------------------------------------------
%   (c) 2025 Technische Universität Wien and Nokia
%   Authored by: Bashar Tahir
%   Project: Christian Doppler Laboratory for Digital Twin assisted AI for sustainable Radio Access Networks
% -----------------------------------------------------------------------
%   Mode-m unfolding of a tensor
%
%%
function Xm = tunfold(X, m)
    sizeArray    = size(X);
    order        = 1:length(sizeArray);
    order(1)     = m;
    order(m)     = 1;
    sizem        = sizeArray(m);
    sizeArray(m) = [];
    Xm = reshape(permute(X, order), [sizem, prod(sizeArray)]);
end

