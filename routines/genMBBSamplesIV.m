function res = genMBBSamplesIV(VAR, UBlocks, ZBlocks, uCenter, zCenter, blkSize, J)
% This function construct moving-block bootstrapped samples for SVAR-IV
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
bootU = zeros(J*blkSize,n);
bootZ = zeros(J*blkSize,1);
for j = 1:J
    bootU(1+blkSize*(j-1):blkSize*j,:) = UBlocks(:,:,index(j,1));
    bootZ(1+blkSize*(j-1):blkSize*j,:) = ZBlocks(:,:,index(j,1));
end
bootU = bootU(1:T,:);
bootZ = bootZ(1:T,:);

%check for the case where all of the proxies equal zero
if sum(bootZ == 0) == T
    res.valid = false;
else
    %center the VAR residuals and proxies
    bootU = bootU - uCenter;
    bootZ(bootZ(:,1)~=0,1) = bootZ(bootZ(:,1)~=0,1) -...
        zCenter(bootZ(:,1)~=0,1);

    for t=1:T

        Ylag =flipud(bootY(t:t+p-1,:))';
        Ylag =[DET(t+p,:)'; Ylag(:)];

        bootY(t+p,:) = Ylag'*coef + bootU(t,:);

    end

    res.Y = bootY;
    res.Z = bootZ;
    res.valid = true;

end