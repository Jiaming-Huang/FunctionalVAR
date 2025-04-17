function [DATASET, modelSpec] = prepareSimFVAR(y,iv,fcn,gridFcn,p,iDET, nFPCMax, identification, irhor, cLevel, nBoot, blkSize, verbose)
DATASET          = struct;
modelSpec        = struct;

DATASET.data     = [y iv];
DATASET.varsName = {'y', 'iv'};
DATASET.gridFcn  = gridFcn; % 1 by nGrid vector
DATASET.fcn      = fcn;
switch identification
    case 'CHOL'
        modelSpec.varsSel = {'y','Func'};
        modelSpec.varsShock = {'y'};
    case 'InternalIV'
        modelSpec.varsSel = {'iv','Func', 'y'};
        modelSpec.varsShock = {'iv'};
        modelSpec.varsUnitNorm = {'y'};
    case 'ExternalIV'
        modelSpec.varsSel = {'y','Func'};
        modelSpec.varsShock = {'y'};
        modelSpec.Instrument = {'iv'};
end
modelSpec.identification = identification;

modelSpec.p = p;
modelSpec.iDET = iDET;
modelSpec.nFPCMax = nFPCMax;
modelSpec.irhor = irhor;
modelSpec.nBoot = nBoot;
modelSpec.cLevel = cLevel;
modelSpec.blkSize = blkSize;
modelSpec.verbose = verbose;
end