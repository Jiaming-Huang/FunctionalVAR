function [y, expanded_varnames] = constructVARMat(data, fpcs, names, MAP)
% Constructs a VAR matrix from aggregate and functional data.
%
% Inputs:
%   - data  : T x n matrix of observed variables
%   - fpcs  : T x k matrix of functional principal components (FPCs)
%   - names : cell array of variable names (strings), including 'Func' for FPCs
%   - MAP   : containers.Map object mapping variable names to column indices in 'data'
%
% Outputs:
%   - y                 : [n x q] matrix containing selected variables and FPCs
%   - expanded_varnames : cell array of names corresponding to the columns of 'y'

m = numel(names);
output_parts = cell(1, m);
expanded_varnames = {};

for i = 1:m
    if strcmp(names{i}, 'Func')
        k = size(fpcs, 2);
        output_parts{i} = fpcs;
        fpc_names = arrayfun(@(j) sprintf('fpc%d', j), 1:k, 'UniformOutput', false);
        expanded_varnames = [expanded_varnames, fpc_names];
    else
        output_parts{i} = data(:, MAP(names{i}));
        expanded_varnames = [expanded_varnames, names(i)];
    end
end

y = [output_parts{:}];
end