function [res] = scatter3d(x, y, z, labels, ms, mkr, filled, convexhull)
%
% 3D-scatter plot of values in x, y and z, colored by labels.
%
% Input arguments :
% =================
% x : values vector for x-axis
% y : values vector for y-axis
% z : values vector for z-axis
% labels : class assignment as a string vector or numerical vector/matrix
% ms : marker size
% mkr : marker type
% filled : whether markers should be filled (1) or not (0) 
% convexhull : whether convex hulls are drawn (1) around classes or not (0)
%
% Output arguments :
% ==================
% res : structure of color and class assignments
% Produces an inline 3D-scatter plot
%
% Usage :
% =======
% ms = 30;
% mkr = 'o';
% filled = 1;
% scatter3d(x, y, z, labels, ms, mkr, filled, convexhull);
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
% 2021.06.01 : uses the matlab plot function, not scatter anymore
%
% =========================================================================

%% Fail-safe section

% Color assignment acording to classes
res = classcol(labels);

if nargin < 7
    convexhull = 1;
end

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
    
    idx = find(res.classvec == i);
    
    if filled == 1
        res.handles{i} = plot3(x(idx), y(idx), z(idx), 'Color', res.colmat(idx(1),:), 'MarkerFaceColor', res.colmat(idx(1),:), 'Marker', mkr{i}, 'MarkerSize', ms, 'LineStyle','None');
        hold on;
    else
        res.handles{i} = plot3(x(idx), y(idx), z(idx), 'Color', res.colmat(idx(1),:), 'Marker', mkr{i}, 'MarkerSize', ms, 'LineStyle','None');
    end
    
    if convexhull == 1
        k = convhulln([x(idx), y(idx) z(idx)]);
        %plot3(x(idx(k)), y(idx(k)), z(idx(k)),':','Color',res.colmat(idx(1),:),'linewidth',2);
        trisurf(k,x(idx), y(idx), z(idx),'FaceColor',res.colmat(idx(1),:),'FaceAlpha',0.4)
    end
end

box on; axis tight;

end
