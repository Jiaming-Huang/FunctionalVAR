% This script specifies parameters for FVAR(1) and computes analytical IRs

%% Basis Functions Setup
q        = 1;
nBasis   = 2*q + 1;
nGrid    = 400;
gridFcn    = linspace(0, 1, nGrid);
basismat = zeros(nGrid, nBasis);

% First basis is constant
basismat(:, 1) = 1;

% Fourier basis functions
for k = 1:q
    basismat(:, 2*k)   = sqrt(2) * cos(2*pi*k*gridFcn);
    basismat(:, 2*k+1) = sqrt(2) * sin(2*pi*k*gridFcn);
end

%% SVAR Parameter Specification
if PARAMS_SET == 1
    A = [ 0.4000  -0.3000   0.2000   0.1000;
         -0.6000  -0.0500  -0.2300   0.7600;
         -0.3000   0.8000  -0.0500   0.0400;
         -0.4000   0.0400   0.7600   0.2300 ];
     
    B = [ 1.0000   0       0       0;
          0.5000   1.0000  0       0;
          0.5000   0.5000  1.0000  0;
          0.5000   0.5000  0.5000  1.0000 ];
      
elseif PARAMS_SET == 2
    A = [ 0.8000  -0.3000   0.2000   0.1000;
          0.2000   0.3000   0.2300   0.7000;
         -0.3000  -0.5000  -0.5000   0.0400;
          0.1000   0.0400  -0.3000   0.6000 ];

    B = [ 1.0000   0       0       0;
          0.3333   1.0000  0       0;
         -0.1000   0       0.5000  0;
          0.0500   0       0       0.3333 ];
end

%% Analytical Impulse Responses
irs0_fpcs = zeros(nBasis, irhor + 1);
irs0_f    = zeros(nGrid, irhor + 1);
for h = 0:irhor
    tmp = A^h * B(:,1);                    % IR of the first shock
    irs0_fpcs(:, h + 1) = tmp(2:end);      % Responses of FPC scores
    irs0_f(:, h + 1)    = basismat * irs0_fpcs(:, h + 1); % Functional IRs
end

%% Plotting

% Plot basis functions
figure;
for j = 1:nBasis
    subplot(q+1, 2, j);
    plot(gridFcn, basismat(:, j));
    xlabel('u');
    ylabel('Basis value');
    title(['Basis ', num2str(j)]);
end
sgtitle('Fourier Basis Functions');

% Plot true functional impulse responses
[gridPlotIrsTime, gridPlotIrsFcn] = meshgrid(0:irhor, gridFcn);

figure;
mesh(gridPlotIrsTime, gridPlotIrsFcn, irs0_f);
xlabel('Horizon');
ylabel('u');
zlabel('Response');
title('True Functional Impulse Responses');