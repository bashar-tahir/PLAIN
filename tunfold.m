%% ------------------------------------------------------------------------
%   (c) 2025 Bashar Tahir, bashar.tahir@tuwien.ac.at
%   Institute of Telecommunications, TU Wien
%   https://www.tuwien.at/etit/tc/en/
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

