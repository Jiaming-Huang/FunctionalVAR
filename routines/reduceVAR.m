function VAR = reduceVAR(data, p, iDET)
% estimate reduced-form VAR 
% ====================================================================%
% Input : 
%   data - T by n data
%   p    - number of lags
%   iDET - indicator for deterministic terms
%           1 for constants, 2 deterministic trend, 3 quadratic trend
% Output : VAR - struct
%   data - raw input, useful for bootstrap
%   Y    - LHS
%   X    - RHS
%   DET  - deterministic trend used
%   p    - number of lags
% ====================================================================%

% prepare variables
T           = size(data,1);
DET         = [ones(T,1) (1:T)' (1:T).^2'];
DET         = DET(:,iDET);

Y           = data;
X           = lag(data,p);

% remove missing values
Y           = Y(p+1:end,:);
X           = [DET(p+1:end,:) X(p+1:end,:)];

[VAR.T, VAR.n] = size(Y);

% OLS
VAR.b       = (X'*X)\(X'*Y);
VAR.u       = Y - X*VAR.b;
VAR.SIGMA   = (VAR.u'*VAR.u)/(VAR.T-size(VAR.b,1)); % B0^{-1}*B0^{-1}'

% outpout
VAR.data    = data;
VAR.Y       = Y;
VAR.X       = X;
VAR.DET     = DET;
VAR.p       = p;
end