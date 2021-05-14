function [res] = elliscatter2d(x, y, labels, alpha, ms, mkr)
%
% 2D-scatter plot of values in x and y, colored by labels.
% Adds confidence ellipse based on reference below.
% 
% Reference :
% ===========
% https://groups.google.com/forum/#!topic/comp.soft-sys.matlab/U8229nQbnts
% (Based on post from Doug Schwarz accessed June 19, 2020)
%
% Input arguments :
% =================
% x : values vector for x-axis 
% y : values vector for y-axis 
% labels : class assignment as a string vector or numerical vector/matrix
% alpha : confidence interval as float between 0 and 1, default is 0.95
% ms : marker size
% mkr : marker type
%
% Output arguments :
% ==================
% res : structure of color and class assignments 
% Produces an inline 2D-scatter plot
%
% Usage :
% =======
% elliscatter2d(x, y, labels, 0.95, 30, 'o');
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
Y = res.classvec; % vector of class assignments

if nargin < 6
    mkr = distmkr(max(Y));
end

if length(mkr) == 1
    mkr = repelem(mkr,max(Y));
elseif length(mkr) ~= max(Y)
    mkr = distmkr(max(Y));
end

if nargin < 5
    ms = 20;
end

if nargin < 4
    alpha = 0.5;
end

if nargin < 3 
    error('At least x and y axes coordinates must be provided, along with class assignment information');
end

if min(size(x)) > 1 || min(size(y)) > 1
    error('Check dimensions of x and y : only vectors are allowed');
end

if iscolumn(x) == 0
    x = x';
end

if iscolumn(y) == 0
    y = y';
end

%% Main section

% Merges two axes
X = [x,y];

[n,p] = size(X); % data dimensions

npoints = 100; % number of points to plot the ellipses
coord = linspace(0,2*pi,npoints)';

Options.prepro.X.type = {'cmeancenter'}; % preprocessing option
Options.pcatype = 'svd'; % pca by svd

F = finv(alpha, p, n-p) * ((p*(n-1))/((n-p))); % confidence interval according to alpha
    
% figure; hold on;

for i = 1 : length(res.labels)
    
    centroid = mean(X(Y == i,:),1); % centroid coordinates
    
    [prepro] = xpreproc(X(Y == i,:), Options);
    [pca_model] = pcac(prepro.data, 2, Options);
    
    evals = pca_model.evals; % eigenvalues
    P = pca_model.loadings; % loadings
    
    ev = diag(sqrt(F*evals));
    exy = [cos(coord),sin(coord)]*ev*P' + repmat(centroid,npoints,1);
    
    res.handles{i} = scatter(x(res.classvec == i), y(res.classvec == i), ms, res.colmat(res.classvec == i,:), 'filled', mkr{i});
    line(exy(:,1),exy(:,2),'color',res.labelscol(i,:),'LineWidth', 2); % adds ellipse to plot
    
end

%legend(res.labels);
box on; axis tight;


end