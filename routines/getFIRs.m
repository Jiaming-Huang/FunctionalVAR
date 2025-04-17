function [irs_agg, irs_f] = getFIRs(VAR, Bzero, bas, isShockVar, isDiffVar, isFcnCols, irhor)
    I           = eye(VAR.n);
    shockMat    = I(:,isShockVar);
    irs_temp    = getIRs(VAR, Bzero, irhor, shockMat, isDiffVar);
    irs_agg     = irs_temp(~isFcnCols,:);
    irs_f       = bas*irs_temp(isFcnCols,:);
end