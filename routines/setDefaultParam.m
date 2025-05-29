function [DATASET, modelSpec] = setDefaultParam(DATASET, modelSpec)
% Set default parameters if missing

% Aggregate data
[T, n]   = size(DATASET.data);
varsName = arrayfun(@(k) sprintf('y%d', k), 1:n, 'UniformOutput', false);
DATASET  = setParam(DATASET, 'varsName', varsName);
DATASET  = setParam(DATASET, 'varsDiff', {});
DATASET  = setParam(DATASET, 'dates', 1:T);

% Functional data
basisobj = create_fourier_basis([0,1], 21);
fdParobj = fdPar(basisobj);
modelSpec = setParam(modelSpec, 'fdParobj', fdParobj);

% VAR
modelSpec = setParam(modelSpec, 'p', 1);
modelSpec = setParam(modelSpec, 'iDET', 1);
modelSpec = setParam(modelSpec, 'dateStart', DATASET.dates(1));
modelSpec = setParam(modelSpec, 'dateEnd', DATASET.dates(end));

% Functional PCA
modelSpec = setParam(modelSpec, 'nFPCMax', 20);

% SVAR & IRs
modelSpec = setParam(modelSpec, 'identification', 'CHOL');
modelSpec = setParam(modelSpec, 'varsShock', modelSpec.varsSel(1));
modelSpec = setParam(modelSpec, 'varsUnitNorm', modelSpec.varsShock);
modelSpec = setParam(modelSpec, 'Instrument', {});
modelSpec = setParam(modelSpec, 'nZlags', 0);
modelSpec = setParam(modelSpec, 'nNWlags', -1);
modelSpec = setParam(modelSpec, 'irhor', 36);
modelSpec.varsDiff  = DATASET.varsDiff;

% Bootstrap & Inference
modelSpec = setParam(modelSpec, 'nBoot', 500);
modelSpec = setParam(modelSpec, 'blkSize', 16);
modelSpec = setParam(modelSpec, 'cLevel', 95);
modelSpec = setParam(modelSpec, 'verbose', false);

end

function s = setParam(s, field, default)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = default;
end
end