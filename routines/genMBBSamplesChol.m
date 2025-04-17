function bootY = genMBBSamplesChol(VAR, UBlocks, uCenter, blkSize, J)
% This function construct moving-block bootstrapped samples with SVAR
% used with caution, determinisitc term is placed first in VAR

% unpack variables
T       = VAR.T;
n       = VAR.n;
p       = VAR.p;
coef    = VAR.b;



%% build bootsrap sample
% initialization
DET     = VAR.DET; % DET is of the same dim as VAR.data!
bootY   = VAR.data;
bootY(p+1:end,:) = nan;


%draw bootstrapped VAR residuals and proxies
index = ceil((T - blkSize + 1)*rand(J,1));
bootUy = zeros(J*blkSize,n);
for j = 1:J
    bootUy(1+blkSize*(j-1):blkSize*j,:) = UBlocks(:,:,index(j,1));
end
bootUy = bootUy(1:T,:);
bootUy = bootUy - uCenter;

for t=1:T

    Ylag =flipud(bootY(t:t+p-1,:))';
    Ylag =[DET(t+p,:)'; Ylag(:)];
    bootY(t+p,:) = Ylag'*coef + bootUy(t,:);

end

end