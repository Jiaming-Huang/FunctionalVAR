function res = getBandsBoot(VAR, FPCAOut, modelSpec)
% calculates bands around IRFs using moving block bootstrap

% unpack variables
% VAR
p           = VAR.p;
T           = VAR.T;
u           = VAR.u;
iDET        = modelSpec.iDET;
varsDiff    = modelSpec.varsDiff;
isFcnCols   = modelSpec.isFcnCols;
nAgg        = length(modelSpec.aggSel);
mapSub      = modelSpec.mapSub;
varsSel     = modelSpec.varsSel;

% FPCA
fdParobj    = modelSpec.fdParobj;
nFPCMax     = modelSpec.nFPCMax;
gridFcn       = FPCAOut.gridFcn;
mu          = FPCAOut.mu;
bas         = FPCAOut.bas;
resFcn      = FPCAOut.resFcn;
ngridFcn      = size(gridFcn,2);
Tfcn        = size(resFcn,1);
% SVAR and IRs
varsShock   = modelSpec.varsShock;
nShocks     = length(varsShock);
irhor       = modelSpec.irhor;

% Bootstrap and Inference
cLevel      = modelSpec.cLevel;
nBoot       = modelSpec.nBoot;
verbose     = modelSpec.verbose;
blkSize     = modelSpec.blkSize;
J           = ceil(T/blkSize);


% dataholder
bootirs     = nan(nAgg,irhor+1,nShocks,nBoot);
bootirs_f   = nan(ngridFcn,irhor+1,nShocks,nBoot);

switch modelSpec.identification
    case 'CHOL'
        %% create the blocks for the moving block bootstrap
        [UBlocks, uCenter] = prepareBlocks(u, blkSize, T, J);

        %the bootstrap iterations
        for iBoot = 1:nBoot
            % construct bootstrap samples
            bootY   = genMBBSamplesChol(VAR, UBlocks, uCenter, blkSize, J);
            % get the bootstrapped functions
            bSwitch = 1-2*(rand(Tfcn,1)>.5);
            bootFcn = mu + bootY(:,isFcnCols)*bas' + resFcn.*bSwitch;

            %% Re-estimation on bootstrapped data
            %% 1 FPCA
            [~, bootbas, bootfpcs]     = doFPCA(gridFcn, bootFcn, fdParobj, nFPCMax);
            % Construct VAR Matrix
            [bootY, expanded_names] = constructVARMat(bootY(:,~isFcnCols), bootfpcs, varsSel, mapSub);
            % dimension of the VAR may have changed in bootstrap!!
            bootIsFcnCols   = contains(expanded_names, 'fpc');
            bootIsShockVar  = ismember(expanded_names,varsShock);
            bootIsDiffVar   = ismember(expanded_names,varsDiff);
            %% 2 SVAR
            bootVAR     = reduceVAR(bootY, p, iDET);
            bootBzero   = chol(bootVAR.SIGMA);
            bootBzero   = bsxfun(@rdivide,bootBzero,diag(bootBzero))';

            %% get IRs
            [bootirs(:,:,:,iBoot), bootirs_f(:,:,:,iBoot)] = getFIRs(bootVAR, bootBzero, bootbas, bootIsShockVar, bootIsDiffVar, bootIsFcnCols, irhor);

            if verbose && (mod(iBoot,100)==0)
                fprintf(1,'bootstrap irf: %i of %i \n', iBoot, nBoot);
            end
        end
        nValidBS = nBoot;
    case 'InternalIV'
        varsUnitNorm        = modelSpec.varsUnitNorm;
        [UBlocks, uCenter]  = prepareBlocks(u, blkSize, T, J);

        %the bootstrap iterations
        for iBoot = 1:nBoot
            % construct bootstrap samples
            bootY   = genMBBSamplesChol(VAR, UBlocks, uCenter, blkSize, J);
            % get the bootstrapped functions
            bSwitch = 1-2*(rand(Tfcn,1)>.5);
            bootFcn = mu + bootY(:,isFcnCols)*bas' + resFcn.*bSwitch;

            %% Re-estimation on bootstrapped data
            %% 1 FPCA
            [~, bootbas, bootfpcs]     = doFPCA(gridFcn, bootFcn, fdParobj, nFPCMax);
            % Construct VAR Matrix
            [bootY, expanded_names] = constructVARMat(bootY(:,~isFcnCols), bootfpcs, varsSel, mapSub);
            % dimension of the VAR may have changed in bootstrap!!
            bootIsFcnCols   = contains(expanded_names, 'fpc');
            bootIsShockVar  = ismember(expanded_names,varsShock);
            bootIsDiffVar   = ismember(expanded_names,varsDiff);

            %% 2 SVAR
            bootVAR         = reduceVAR(bootY, p, iDET);
            bootBzero       = chol(bootVAR.SIGMA);
            bootBzero       = bsxfun(@rdivide,bootBzero,diag(bootBzero))';
            bootIsUnitNorm  = ismember(expanded_names, varsUnitNorm);
            bootBzero(:,bootIsShockVar) = bootBzero(:,bootIsShockVar)/bootBzero(bootIsUnitNorm,bootIsShockVar);


            %% get IRs
            [bootirs(:,:,:,iBoot), bootirs_f(:,:,:,iBoot)] = getFIRs(bootVAR, bootBzero, bootbas, bootIsShockVar, bootIsDiffVar, bootIsFcnCols, irhor);

            if verbose && (mod(iBoot,100)==0)
                fprintf(1,'bootstrap irf: %i of %i \n', iBoot, nBoot);
            end
        end
        nValidBS = nBoot;

    case 'ExternalIV'
        % some params specific to external IV
        nZlags      = modelSpec.nZlags;
        nNWlags     = modelSpec.nNWlags;
        validBS     = true(nBoot,1);
        %% create the blocks for the moving block bootstrap
        [UBlocks, ZBlocks, uCenter, zCenter] = prepareBlocksIV(u, modelSpec.z, blkSize, T, J);

        for iBoot = 1:nBoot
            % construct bootstrap samples
            bootSet = genMBBSamplesIV(VAR, UBlocks, ZBlocks, uCenter, zCenter, blkSize, J);
            if bootSet.valid
                bootY   = bootSet.Y;
                bootZ   = bootSet.Z;
                % get the bootstrapped functions
                bSwitch = 1-2*(rand(Tfcn,1)>.5);
                bootFcn = mu + bootY(:,isFcnCols)*bas' + resFcn.*bSwitch;

                %% Re-estimation on bootstrapped data
                %% 1 FPCA
                [~, bootbas, bootfpcs]     = doFPCA(gridFcn, bootFcn, fdParobj, nFPCMax);
                % Construct VAR Matrix
                [bootY, expanded_names] = constructVARMat(bootY(:,~isFcnCols), bootfpcs, varsSel, mapSub);
                % dimension of the VAR may have changed in bootstrap!!
                bootIsFcnCols   = contains(expanded_names, 'fpc');
                bootIsShockVar  = ismember(expanded_names,varsShock);
                bootIsDiffVar   = ismember(expanded_names,varsDiff);

                %% 2 SVAR
                bootVAR         = reduceVAR(bootY, p, iDET);
                [bootBzero, ~]   = getBzeroIV(bootVAR, bootZ, bootIsShockVar, nZlags, nNWlags, false);
                % bootSmplIV  = ~any(isnan(bootZ),2);
                % bootRes     = ProxySVARidentification(bootVAR.u(bootSmplIV,:), bootIsShockVar, bootZ(bootSmplIV,:),size(bootVAR.b,1));
                % bootBzero   = zeros(bootVAR.n,bootVAR.n);
                % bootBzero(:,bootIsShockVar) = bootRes.B/bootRes.B(bootIsShockVar);

                %% get IRs
                [bootirs(:,:,:,iBoot), bootirs_f(:,:,:,iBoot)] = getFIRs(bootVAR, bootBzero, bootbas, bootIsShockVar, bootIsDiffVar, bootIsFcnCols, irhor);

            else
                % in the very unlucky case, we sample an IV with all
                % zeros, or all nans...
                validBS(iBoot)=false;
            end

            if verbose && (mod(iBoot,100)==0)
                fprintf(1,'bootstrap irf: %i of %i \n', iBoot, nBoot);
            end
        end

        % sanity check
        if sum(validBS)/nBoot <0.5
            warning('Less than half valid Bootstrap samples!');
        end
        bootirs      = bootirs(:,:,:,validBS);
        bootirs_f    = bootirs_f(:,:,:,validBS);
        nValidBS     = sum(validBS);
end


% get relevant quantiles
uBound = cLevel+(100-cLevel)/2;
uBound = uBound/100;
lBound = (100-cLevel)/2;
lBound = lBound/100;

% irs for aggregates
bootirs_f = sort(bootirs_f,4);
res.upper_f = bootirs_f(:,:,:,round(uBound*nValidBS));
res.lower_f = bootirs_f(:,:,:,round(lBound*nValidBS));

% functional irs
bootirs =sort(bootirs,4);
res.upper_agg =bootirs(:,:,:,round(uBound*nValidBS));
res.lower_agg =bootirs(:,:,:,round(lBound*nValidBS));
res.nValidBS = nValidBS;
end


function [blocks, center] = prepareBlocks(data, blkSize, T, J)
    nBlks       = T-blkSize+1;
    blocks = zeros(blkSize, size(data,2), nBlks);
    for ii = 1:nBlks
        blocks(:,:,ii) = data(ii:blkSize+ii-1,:);
    end
    %center the bootstrapped VAR residuals and the IV
    center = mean(blocks,3);
    center = repmat(center,[J,1]);
    center = center(1:T,:);
end

function [UBlocks, ZBlocks, uCenter, zCenter] = prepareBlocksIV(u, z, blkSize, T, J)
    [UBlocks, uCenter] = prepareBlocks(u, blkSize, T, J);
    
    % Handle Z blocks with missing values / zeros
    nIV         = size(z,2);
    ZBlocks = zeros(blkSize, nIV, T-blkSize+1);
    for ii = 1:T-blkSize+1
        ZBlocks(:,:,ii) = z(ii:blkSize+ii-1,:);
    end
    
    % Z may contain many zeros and/or NaN, calculate mean excluding
    % those observations
    validMask = (ZBlocks ~= 0) & ~isnan(ZBlocks);
    zClean = ZBlocks;
    zClean(~validMask) = 0;
    zSum = sum(zClean, 3);
    zCount = sum(validMask, 3);
    % Avoid division by zero
    zCount(zCount == 0) = NaN;
    zCenter = zSum ./ zCount;
    zCenter = repmat(zCenter,[J,1]);
    zCenter = zCenter(1:T,:);
end

