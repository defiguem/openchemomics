function [res] = scatter2d(x, y, labels, ms, mkr, filled)
%
% 2D-scatter plot of values in x and y, colored by labels.
% 
% Input arguments :
% =================
% x : values vector for x-axis 
% y : values vector for y-axis 
% labels : class assignment as a string vector or numerical vector/matrix
% ms : marker size
% mkr : marker type
% filled : whether markers should be filled (1) or not (0) 
%
% Output arguments :
% ==================
% res : structure of color and class assignments 
% Produces an inline 2D-scatter plot
%
% Usage :
% =======
% ms = 30;
% mkr = 'o';
% filled = 1;
% scatter2d(x, y, labels, ms, mkr, filled);
%
% Related functions :
% ===================
% classcol.m (distinguishable colors and class assignment)
% distmkr.m (distinguishable markers)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
%
% =========================================================================

%% Fail-safe section

% Color assignment acording to classes
res = classcol(labels);

if nargin < 6
    filled = 1;
end

if nargin < 5
    mkr = distmkr(max(res.classvec));
end

if ischar(mkr)
    mkr = {mkr};
end

if length(mkr) == 1
    mkr = repelem(mkr,max(res.classvec));
elseif length(mkr) ~= max(res.classvec) || isempty(mkr) == 0
    mkr = distmkr(max(res.classvec));
end

if nargin < 4
    ms = 30;
end

if nargin < 3 
    error('At least x and y axes coordinates must be provided, along with class assignment information');
end

if min(size(x)) > 1 || min(size(y)) > 1
    error('Check dimensions of x and y : only vectors are allowed');
end

%% Main section

for i = 1 : length(res.labels)
    if filled == 1
        res.handles{i} = scatter(x(res.classvec == i), y(res.classvec == i), ms, res.colmat(res.classvec == i,:), 'filled', mkr{i});
        hold on;
    else
        res.handles{i} = scatter(x(res.classvec == i), y(res.classvec == i), ms, res.colmat(res.classvec == i,:), mkr{i});
    end
end

box on; axis tight;

end