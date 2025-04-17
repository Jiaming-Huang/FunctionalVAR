% Simulation: FVAR(1)
%% Housekeeping
clear;
close all;
clc

addpath('../fdaM')
addpath('../routines')

rng(27);

%% 1. Setup
% set params for
% dgp
PARAMS_SET          = 1;
gridT               = [100,200,500];
% estimation
nFPCMax             = 5; % a conservative number of FPC
p                   = 1;
iDET                = 1;
identification      = 'ExternalIV';
irhor               = 12;
nBoot               = 500;
blkSize             = 16;
cLevel              = 95;
verbose             = false;
% simulation
nRep                = 1000;

% generate params
simulSetParamsFVAR;

% output containers
nGridT              = length(gridT);
SIM_IR_FVAR         = nan(nGrid, irhor+1, nRep, nGridT);
SIM_BIAS_FVAR       = nan(nGrid, irhor+1, nRep, nGridT);
SIM_ERRORS_FVAR     = nan(irhor+1, 2, nRep, nGridT);
SIM_UBANDS_FVAR     = nan(nGrid, irhor+1, nRep, nGridT);
SIM_LBANDS_FVAR     = nan(nGrid, irhor+1, nRep, nGridT);
if strcmp(identification,'InternalIV')
    nAgg = 2;
else
    nAgg = 1;
end
SIM_IR_AGG          = nan(nAgg, irhor+1, nRep, nGridT);
SIM_UBANDS_AGG      = nan(nAgg, irhor+1, nRep, nGridT);
SIM_LBANDS_AGG      = nan(nAgg, irhor+1, nRep, nGridT);
irs0_y = zeros(size(A,1), irhor+1);
for h = 0:irhor
    irs0_y(:,h+1) = A^h*B(:,1);
end
irs0_agg = irs0_y(1,:);

%% 2. Simulation 
startGrid = tic;
for tt = 1:nGridT
    T = gridT(tt);
    fprintf('........Working with sample size T= %d...........\n', T);
    for iRep = 1:nRep
        %% Step 0: generate data from VAR(1) with FPC scores
        [Y, iv] = genSVARIV(A, B, T);
        y     = Y(:,1); % aggregate variable
        fpcs0 = Y(:,2:end); % FPC scores
        fcn     = fpcs0*basismat'; % functional data

        % prepare input for doFVAR
        [DATASET, modelSpec] = prepareSimFVAR(y,iv,fcn,gridFcn,p,iDET, nFPCMax, identification, irhor, cLevel, nBoot, blkSize, verbose);
        
        FVAR = doFVAR(DATASET, modelSpec);

        %% Store output
        SIM_IR_FVAR(:,:,iRep,tt) = FVAR.irs_f;
        SIM_UBANDS_FVAR(:,:,iRep,tt) = FVAR.upper_f;
        SIM_LBANDS_FVAR(:,:,iRep,tt) = FVAR.lower_f;
        SIM_IR_AGG(:,:,iRep,tt) = FVAR.irs_agg;
        SIM_UBANDS_AGG(:,:,iRep,tt) = FVAR.upper_agg;
        SIM_LBANDS_AGG(:,:,iRep,tt) = FVAR.lower_agg;

        % show progress
        fprintf('........Replication iRep = %d...........\n', iRep);
        % if mod(iRep/50,1)==0
        %     fprintf('Iteration %d \n', iRep)
        %     fprintf('Mean L2 errors so far: %f \n',mean(mean(SIM_IR_ERRORS(:,1,1:iRep))));
        %     fprintf('Mean uniform errors so far: %f\n',mean(mean(SIM_IR_ERRORS(:,2,1:iRep))));
        % end

    end
end
endGrid = toc(startGrid);
fprintf('Grid finished. Time used: %f seconds.\n', endGrid)

outname = strcat(['../output/simul/fvar_param',num2str(PARAMS_SET),'.mat']);
save(outname)


%% Evaluate FVAR - Aggregate IRs
bias = squeeze(mean(irs0_agg - SIM_IR_AGG,3));
mse = squeeze(mean((irs0_agg - SIM_IR_AGG).^2,3));
figure;
subplot(1,2,1);
plot(0:irhor, bias,'LineWidth',2);
legend(compose('T = %d', gridT));
title('Bias');
subplot(1,2,2);
plot(0:irhor, mse,'LineWidth',2);
legend(compose('T = %d', gridT));
title('MSE');
sgtitle('Aggregate IRs');

SIM_COVERAGE_AGG = squeeze(mean((SIM_UBANDS_AGG >= irs0_agg) & (SIM_LBANDS_AGG <= irs0_agg),3));
figure;
plot(0:irhor, SIM_COVERAGE_AGG,'LineWidth',2);
legend(compose('T = %d', gridT));
title('Coverage');

%% Evaluate FVAR - Functional IRs
[gridPlotIrsTime, gridPlotIrsFcn] = meshgrid(0:irhor,gridFcn);

% bias
bias = squeeze(mean(irs0_f - SIM_IR_FVAR,3));
figure;
for tt = 1:nGridT
    subplot(2,ceil(nGridT/2),tt);
    mesh(gridPlotIrsTime,gridPlotIrsFcn,bias(:,:,tt));
    title(strcat(['T = ',num2str(gridT(tt))]));
end
sgtitle('Bias');

% errors
errs = nan(irhor+1,3,nRep,nGridT);
for tt = 1:nGridT
    for iRep = 1:nRep
        errs(:,1,iRep,tt) = max(abs(irs0_f-SIM_IR_FVAR(:,:,iRep,tt)));
        errs(:,2,iRep,tt) = trapz(gridFcn,abs(irs0_f-SIM_IR_FVAR(:,:,iRep,tt)));
        errs(:,3,iRep,tt) = trapz(gridFcn,(irs0_f-SIM_IR_FVAR(:,:,iRep,tt)).^2);
    end
end
errs = squeeze(mean(errs,3));
figure;
for tt = 1:nGridT
    subplot(2,2,tt);
    plot(0:irhor,errs(:,:,tt),'LineWidth',1.5);
    grid on;
    xlim([0,irhor])
    title(strcat(['T = ',num2str(gridT(tt))]))
    if tt == 1
        legend('UE','IAE','ISE');
    end
end
sgtitle('Estimation Errors');

% coverage
SIM_COVERAGE_FVAR = squeeze(mean((SIM_UBANDS_FVAR >= irs0_f) & (SIM_LBANDS_FVAR <= irs0_f),3));
figure;
for tt = 1:nGridT
    subplot(2,2,tt);
    mesh(gridPlotIrsTime,gridPlotIrsFcn,SIM_COVERAGE_FVAR(:,:,tt));
    title(strcat(['T = ',num2str(gridT(tt))]))
end

coverage_fvar = squeeze(mean(SIM_COVERAGE_FVAR));

% get table
for tt = 1:nGridT
    disp(strcat('T=',num2str(gridT(tt))));
    disp([errs(:,:,tt) coverage_fvar(:,tt)]);
end
