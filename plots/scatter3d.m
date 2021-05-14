function [res] = scatter3d(x, y, z, labels, ms, mkr, filled, type)
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
% scatter3d(x, y, labels, ms, mkr, filled);
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

if nargin < 7
    type = 'qual'; % 'qualitative' as opposed to 'sequential' ('seq')
end

if nargin < 6
    filled = 1;
end

if nargin < 5
    mkr = 'o';
end

if nargin < 4
    ms = 20;
end

if nargin < 3
    error('At least x, y and z axes coordinates muste be provided, along with class assignment information');
end

if min(size(x)) > 1 || min(size(y)) > 1 || min(size(z)) > 1
    error('Check dimensions of x and y : only vectors are allowed');
end

%% Main section

% Color assignment acording to classes
res = classcol(labels);

if strcmp(type,'qual') == 0
    % Defines a gradient of colors between col1 and col2
    col1 = [0, 0.4470,0.7410];
    col2 = [0.9769 0.9839 0.0805];
    [~,I] = sort(labels,'ascend');
    lgth = length(labels);
    grad = [linspace(col1(1),col2(1),lgth)',linspace(col1(2),col2(2),lgth)',linspace(col1(3),col2(3),lgth)'];
    for i = 1 : lgth
        res.colmat(i,:) = grad((I==i),:);
    end
end

figure; hold on;
for i = 1 : length(res.labels)
    if filled == 1
        scatter3(x(res.classvec == i), y(res.classvec == i), ...
            z(res.classvec == i), ms, res.colmat(res.classvec == i,:), 'filled', mkr);
    else
        scatter3(x(res.classvec == i), y(res.classvec == i), ...
            z(res.classvec == i), ms, res.colmat(res.classvec == i,:), mkr);
    end
end

box on; axis tight;
if strcmp(type,'qual') == 0
    colormap(grad);
    cbr = colorbar('northoutside');
    set(cbr,'TickLabels',linspace(floor(min(labels)), round(max(labels)), 11));
end

end
