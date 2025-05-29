function FVAR = doFVAR(DATASET, modelSpec)
% This function estimate a Funtional VAR with moving-block bootstrap.
% It allows for 3 identification schemes: short-run restriction, internal
% IV, and external IV
%
% --------------------------- INPUT --------------------------------
%   DATASET   - struct with full data
%       data       - T x N full data matrix of macro aggregates; only
%                       selected variables (columns) are used in VAR
%                       estimation, specified by
%       varsName   - 1 x N cell array of variable names of data
%       varsDiff   - a cell array of variable names which have been
%                       differenced; for those variables, cumulative IRs
%                       will be computed. Example: empty cell varsDiff =
%                       {}, varsDiff = {CPI}
%       dates      - T x 1 datetime
%       fcn        - T x ngridFcn matrix of functional data, where ngridFcn
%                       is the number of grid points
%
%   modelSpec - struct with model specification
%       ---------------  VAR ----------------
%       varsSel:   - 1 x nSel cell array of variable names to be included
%                       in the VAR estimation; if identification = 'CHOL',
%                       the order will be used for the Cholesky decomposition
%       p:         - the lag order of the VAR
%       iDET:      - an integer or vector of integers indicating the
%                       deterministic terms to be included in the VAR;
%                       1: constant, 2: linear trend, 3: quadratic trend.
%                       Example: modelSpec.iDET = 1;
%                       modelSpec.iDET = [1, 2];
%       dateStart: - sample starting period, can be double (e.g. 1991) or
%                       datetime object
%       dateEnd:   - sample ending period, can be double (e.g. 2007) or
%                       datetime object
%       ---------------  FPCA ----------------
%       fdParobj:  - a fdPar object in the fdaM package, used to transform
%                       discretized functional data into functional data object
%       nFPCMax:   - the maximal number of functional principal components;
%                       the actual FPC used in SVAR will be selected by
%                       variation explained and is smaller than nFPCMax
%       ---------------  SVAR & IRs ----------------
%       identification: - a string, either "CHOL", "InternalIV", or "ExternalIV"
%       varsShock: - a cell array of variable names whose corresponding
%                       structural shocks will be used to compute IRs
%       varsUnitNorm: - a cell array of variable names whose contemporanous response
%                       will be normalized to 1; this is only used for the
%                       "InternalIV" identification scheme
%       Instrument: - a cell array of variable names used as instruments;
%                       this is only used for the "ExternalIV" identification scheme
%       irhor     - number of horizons (we compute h=0 by default), so
%                       actually it's irhor+1 horizons
%       ---------------  Bootstrap & Inference ----------------
%       nBoot:     - the number of bootstrap replications
%       blkSize:   - the size of the moving block for bootstrap
%       cLevel:    - the confidence level for the bootstrap inference;
%                       by default, it is 0.95
%       verbose:   - a logical value indicating whether to print the
%                       progress of the bootstrap; by default, it is false
%
% --------------------------- OUTPUT --------------------------------
%   FVAR      - struct with with IRs, bands, and diagnostics

%% unpack params
[DATASET, modelSpec] = setDefaultParam(DATASET, modelSpec);

% Sample
dateStart = modelSpec.dateStart;
dateEnd   = modelSpec.dateEnd;

% Functional PCA
fdParobj  = modelSpec.fdParobj;
nFPCMax   = modelSpec.nFPCMax;

% VAR
varsSel   = modelSpec.varsSel; % the ordering matters (with short-run identification)
varsShock = modelSpec.varsShock;
p         = modelSpec.p;
iDET      = modelSpec.iDET;
identification = modelSpec.identification;
irhor     = modelSpec.irhor;

% extract sampling periods
dates     = DATASET.dates;
ismpl     = dates >= dateStart & dates <= dateEnd;
dateVAR   = dates(ismpl);

% create a map to extract variables later
mapFull   = containers.Map(DATASET.varsName, 1:numel(DATASET.varsName));

%% Step 1: Functional PCA
gridFcn     = DATASET.gridFcn; % 1 by nGridFcn vector
fcn         = DATASET.fcn; % T by nGrid matrix
[mu, bas, fpcs] = doFPCA(gridFcn, fcn, fdParobj, nFPCMax);
fhat        = mu + fpcs * bas';
resFcn      = fcn - fhat; % FPCA residuals
FPCAOut     = struct('gridFcn', gridFcn, 'bas', bas, 'mu', mu, 'resFcn', resFcn);

%% Step 2: SVAR
%% 2.1 Prepare data matrix for VAR
% extract aggregate variables
data                = DATASET.data(ismpl, :);
[y, expanded_names] = constructVARMat(data, fpcs, varsSel, mapFull);

% for later usage, redefine the mapping
modelSpec.aggSel = varsSel(~strcmp(varsSel,'Func'));
modelSpec.mapSub = containers.Map(modelSpec.aggSel,1:numel(modelSpec.aggSel));

% check if there are any missing values, if so, modify sampling periods
naDetect = any(isnan(y),2);
if any(naDetect)
    y       = y(~naDetect,:);
    dateVAR = dateVAR(~naDetect);
    fprintf('Your data contains missing values.\nTime range updated to [%s, %s].\n', ...
        datestr(dateVAR(1)), datestr(dateVAR(end)));
end

%% 2.2 Estimate Reduced-form VAR (Cholesky)
VAR = reduceVAR(y, p, iDET);

%% 2.3 Identification
isShockVar = ismember(expanded_names, varsShock);

switch identification
    case 'CHOL'
        %% SHORT-RUN IDENTIFICATION
        % ut = B0^{-1}*wt ==>  Sig_u = B0^{-1}*B0^{-1}'
        % SR restrictions impose that B0 (and thus B0^{-1}) is LOWER-TRIANGULAR
        % So, ignoring the lags, we have
        % Y1t = b11 * ep1t
        % Y2t = b21 * ep1t + b22 * ep2t
        % ...
        % YKt = bK1 * ep1t + ... + bKK * epKt

        Bzero = chol(VAR.SIGMA);   % actually, it's B0^-1'
        Bzero = bsxfun(@rdivide,Bzero,diag(Bzero))'; % remember to transpose Bzero here

    case 'InternalIV'
        % same as Cholesky, except that usually we need normalization for
        % interpretation
        isUnitNorm = ismember(expanded_names, modelSpec.varsUnitNorm);
        Bzero      = chol(VAR.SIGMA);   % actually, it's B0^-1'
        Bzero      = bsxfun(@rdivide,Bzero,diag(Bzero))'; % remember to transpose Bzero here
        Bzero(:,isShockVar) = Bzero(:,isShockVar)/Bzero(isUnitNorm,isShockVar);

    case 'ExternalIV'
        % external IV specific params
        nZlags      = modelSpec.nZlags;
        nNWlags     = modelSpec.nNWlags;
        % extract IV from the dataset
        dateIVStart = dateVAR(p+1);
        ismpl       = dates >= dateIVStart & dates <= dateVAR(end);
        z           = DATASET.data(ismpl, cell2mat(values(mapFull, modelSpec.Instrument)));
        % run 2SLS
        [Bzero, Fstat]  = getBzeroIV(VAR, z, isShockVar, nZlags, nNWlags, true);
        % save for later usage in bootstrap
        modelSpec.z       = z;

    otherwise
        error('Unknown identification method.');
end

%% Step 3: Recover Functional IRs
% construct shockMat given index of the shock
isDiffVar       = ismember(expanded_names, modelSpec.varsDiff);
isFcnCols       = contains(expanded_names, 'fpc');
[irs_agg, irs_f]= getFIRs(VAR, Bzero, bas, isShockVar, isDiffVar, isFcnCols, irhor);

%% Compute Confidence Bands
modelSpec.isFcnCols = isFcnCols;
bootOut = getBandsBoot(VAR, FPCAOut, modelSpec);

%% Output Struct
FVAR            = struct;
FVAR.irs_agg    = irs_agg;
FVAR.upper_agg  = bootOut.upper_agg;
FVAR.lower_agg  = bootOut.lower_agg;
FVAR.irs_f      = irs_f;
FVAR.upper_f    = bootOut.upper_f;
FVAR.lower_f    = bootOut.lower_f;

FVAR.reduceForm = VAR;
FVAR.Bzero      = Bzero;
FVAR.modelSpec  = modelSpec;

if strcmp(identification, 'ExternalIV')
    FVAR.Fstat    = Fstat;
    FVAR.nValidBS = bootOut.nValidBS;
end

end

