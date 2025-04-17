function irs = getIRs(VAR, Bzero, irhor, shockMat, iDiff)
% This function compute impulse responses given reduced form VAR and B0^{-1}'
% Reference: Kilian and Lutkepohl (2017) p110
% 
% --------------------------- INPUT --------------------------------
%   VAR        - struct, output from reduceVAR 
%   Bzero      - B0^{-1}
%   irhor      - integer, impulse response horizons
%   shockMat   - n x nShock matrix, each column consists of the shock
%                   vector
%   iDiff      - vector, whether the VAR variables have been differenced,
%                   if so, compute cumulative IRs
% --------------------------- OUTPUT --------------------------------
%   irs        - n by H by nShock IRFs


%% STEP 0: UNPACK VARIABLES
n           = VAR.n;
p           = VAR.p;
coef        = VAR.b;
nShock      = size(shockMat,2);

%% STEP 1: CONSTRUCT COMPANION FORM (eq 4.1.1)
A           = zeros(n*p, n*p);
A(1:n,:)    = coef(size(VAR.DET,2)+1:end,:)'; 
A(n+1:end, 1:n*(p-1)) = eye(n*(p-1));

%% STEP 2: COMPUTE IRFs
irs      = nan(n,irhor+1,nShock);
for jj = 1:nShock
    irs(:,1,jj) = Bzero*shockMat(:,jj);   
    
    for h = 1:irhor
        Ah = A^h;
        irs(:,h+1,jj) = Ah(1:n,1:n)*Bzero*shockMat(:,jj);   
    end
    
    % obtain cumulative responses if the variables are first-differenced
    for i = 1: n
        if iDiff(i) == 1
            irs(i,:,jj) = cumsum(irs(i,:,jj));
        end
    end
end

end