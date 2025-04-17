function [mu, bas, fpcs] = doFPCA(gridFcn, fcn, fdParobj, nFPCMax)
    %% 1 FPCA
    fdobj = smooth_basis(gridFcn',  fcn', fdParobj);
    pcaOut  = pca_fd(fdobj, nFPCMax, fdParobj);
    
    % extract the output
    mu      = eval_fd(gridFcn, pcaOut.meanfd)';
    bas     = eval_fd(gridFcn, pcaOut.harmfd); % eigenfunctions
    fpcs    = pcaOut.harmscr; % FPC scores
    
    % select the number of FPC based on variation explained >= 95%
    nFPC    = find(cumsum(pcaOut.varprop)>=0.95,1);
    bas     = bas(:,1:nFPC);
    fpcs    = fpcs(:,1:nFPC);
end