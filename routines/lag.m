function [xlag] = lag(x,n)
%take lags/leads of a variable x (t by K)
%   n must be an integer
assert(mod(n, 1) == 0, '%s: n must be an integer', 'lag.m');
xlag = [];
if n>=0
    for i=1:n
        tmp = nan(size(x));
        tmp(i+1:end,:) =  x(1:end-i,:);
        xlag = [xlag tmp];
    end
else
    for i=-1:-1:n
        tmp = nan(size(x));
        tmp(1:end+i,:) = x(1-i:end,:);
        xlag = [xlag tmp];
    end
end

end