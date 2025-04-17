T  = 100;
PARAMS_SET  = 1;
simulSetParamsFVAR;
[Y, iv] = genSVARIV(A, B, T);
y       = Y(:,1); % aggregate variable
fpcs0   = Y(:,2:end); % FPC scores
fcn     = fpcs0*basismat'; % functional data

% prepare input for doFVAR
DATASET = struct;
DATASET.data = [y iv];
DATASET.gridFcn = gridFcn;
DATASET.fcn = fcn;
DATASET.varsName = {'y', 'iv'};
modelSpec = struct;
modelSpec.varsSel = {'iv','Func','y'};
modelSpec.irhor = irhor;
modelSpec.identification = 'InternalIV';
modelSpec.varsUnitNorm = {'y'};

% Implement FVAR
FVAR = doFVAR(DATASET, modelSpec);

% Plotting
[gridPlotIrsTime, gridPlotIrsFcn] = meshgrid(0:irhor,gridFcn);
figure;
subplot(1,2,1);
mesh(gridPlotIrsTime,gridPlotIrsFcn,irs0_f);
title('IR True (Design 1)')
subplot(1,2,2);
mesh(gridPlotIrsTime,gridPlotIrsFcn,FVAR.irs_f);
title('IR Estimated')
figure;
for i=1:irhor
    subplot(3,4,i)
    plot(gridFcn,irs0_f(:,i),'k-','LineWidth',1.5);
    hold on;
    plot(gridFcn,FVAR.irs_f(:,i),'b-');
    plot(gridFcn,FVAR.upper_f(:,i),'r-');
    plot(gridFcn,FVAR.lower_f(:,i),'r-');
end
legend('True','Estimated');