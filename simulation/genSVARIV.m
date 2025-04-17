function [Y, iv] = genSVARIV(A, B, T)
% This function simulates data from a VAR(1) model with an instrument
%   Y_t = A * Y_{t-1} + B * e_t
%   z_t = e_{1t} + v_t
%   e_{1t}, v_t are iid N(0,1)

nBurn = 200;
% initialization
n    = size(A,1);
Y    = zeros(n, T+nBurn);
ep   = randn(n, T+nBurn);

for t = 2:T+nBurn
    Y(:,t) = A * Y(:,t-1) + B * ep(:,t);
end
Y  = Y(:,nBurn+1:end)';
ep = ep(:,nBurn+1:end)';

iv = ep(:,1) + randn(T,1);

end