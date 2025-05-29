function [Bzero, Fstat] = getBzeroIV(VAR, z, isShockVar, nZlags, nNWlags, getFstat)
%GETBZEROIV Computes structural impact vector B and first-stage F-statistic in SVAR-IV
%
%   This function estimates the impact vector in SVAR-IV models where a
%   structural shock is instrumented using an external instrument.
%
%   Reference: Stock and Watson (2018), *Journal of Economic Perspectives*, pp. 1215-1219.
%
%   INPUT:
%       VAR         - structure with fields:
%                       Y: (T x n) matrix of VAR endogenous variables
%                       X: (T x k) matrix of VAR RHS regressors
%                       n: scalar, number of endogenous variables
%       z           - (T x 1) instrument time series
%       isShockVar  - scalar, index of the endogenous variable receiving the structural shock
%       nZlags      - number of lags of the instrument z to include as controls
%       nNWlags     - (optional) number of Newey-West lags; if empty, set automatically
%       getFstat    - logical, whether compute first-stage F statistics
%                       with finite-sample adjustment, this is the same as
%                       running ivreg2 pol Z, robust bw(nNWlags) small and
%                       test the coefficient of z
%
%   OUTPUT:
%       Bzero       - (n x n) impact matrix such that u_t = Bzero * eps_t
%       Fstat       - first-stage F-statistic for instrument relevance
%
% -------------------------------------------------------------------------

% Extract variables
Y = VAR.Y;
pol = VAR.Y(:, isShockVar);               % Shocked policy variable
w = [VAR.X, lag(z, nZlags)];              % Controls: VAR regressors + lags of instrument

% Drop rows with any NaN values (to align observations)
validIdx = ~any(isnan([Y, pol, z, w]), 2);
Y = Y(validIdx, :);
pol = pol(validIdx);
z = z(validIdx);
w = w(validIdx, :);

% Define regressor and instrument matrices
X = [pol, w];                             % Endogenous regressor matrix
Z = [z, w];                               % Instrument matrix

% First stage
Pi = (Z' * Z) \ (Z' * X);                 % First-stage projection coefficients

% Second stage: instrumented regression
Xhat = Z * Pi;                            % Fitted regressor
bhat = (Xhat' * Xhat) \ (Xhat' * Y);      % IV estimates

% Construct B matrix: only one column is identified
Bzero = zeros(VAR.n);
Bzero(:, isShockVar) = bhat(1, :)';


if getFstat
    L = size(z,2);
    [T, P] = size(Z);
    % Set automatic Newey-West lags if not specified
    if nNWlags == -1
        nNWlags = max(0, round(1.3 * sqrt(T)));
    end
    uhat = X(:,1) - Z * Pi(:,1);              % First-stage residuals
    Zu = Z .* uhat;     % For HAC estimation
    v = HACNW(Zu, nNWlags);                      % HAC estimator of covariance
    vPi = (Z' * Z) \ v / (Z' * Z) * T/(T-P);          % HAC variance of Pi, with finite-sample adjustment
    R = zeros(L,P);
    R(1:L,1:L) = eye(L);
    RPi = R * Pi(:,1);
    Fstat = RPi'*((R*vPi*R')\RPi) / L;
else
    Fstat = nan;
end

end