%% ------------------------------------------------------------------------
%   (c) 2025 Bashar Tahir, bashar.tahir@tuwien.ac.at
%   Institute of Telecommunications, TU Wien
%   https://www.tuwien.at/etit/tc/en/
% -----------------------------------------------------------------------
%   Mode-m tensor multiplication.
%   tmul(X, {U_1, U_2, ..., U_M}) calculates mode-m products of X with the
%   factors U_1, U_2, ..., U_M, respectively. Additionally, a certain
%   operation can be applied to the factors before multiplication. See below.
%
%%
function X = tmul(X, U, op)
    if ~exist("op","var") || isempty(op)
        order = 1:length(U);
        for m = 1:length(U)
            orderin = order;
            orderin(1) = m;
            orderin(m) = 1;
            X = permute(pagemtimes(U{m}, permute(X, orderin)), orderin);
        end
    else
        switch op
            case 'T'
                order = 1:length(U);
                for m = 1:length(U)
                    orderin = order;
                    orderin(1) = m;
                    orderin(m) = 1;
                    X = permute(pagemtimes(U{m}.', permute(X, orderin)), orderin);
                end
            case 'H'
                order = 1:length(U);
                for m = 1:length(U)
                    orderin = order;
                    orderin(1) = m;
                    orderin(m) = 1;
                    X = permute(pagemtimes(U{m}', permute(X, orderin)), orderin);
                end
            case 'pinv'
                order = 1:length(U);
                for m = 1:length(U)
                    orderin = order;
                    orderin(1) = m;
                    orderin(m) = 1;
                    X = permute(pagemtimes(pinv(U{m}), permute(X, orderin)), orderin);
                end
        end
    end
end

