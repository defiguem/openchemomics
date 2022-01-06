function [avgmat, anovamat, ssq, ssqvarexp] = sadecomp(X, Ypd, adecomp_model, Options)
%
% The function sadecomp stands for Synthetic ANOVA decomposition. It
% performs the multivariate ANalysis Of VAriance decomposition of the X 
% matrix according to a previous ANOVA decomposition as reference. 
%
% This function was first designed to be used when a resampling procedure
% is used before performing an ANOVA decomposition. It produces the mean
% effect matrices according to the ANOVA model along with the ANOVA
% matrices.
%
% Reference :
% ===========
%
% Input arguments :
% =================
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Ypd : matrix of factors in the columns and levels for each sample
%     identified by integers (n x q), the prediction design matrix
%
% adecomp_model : the ANOVA "model" built using the adecomp.m function
%
% Output arguments :
% ==================
% avgmat : array of cells containing mean factor/interaction effects 
% anovamat : array of cells containing ANOVA-deflated matrices 
%       contained in cells 
% ssq : sum of squares for each effect based on avgmat
% ssqvarexp :  variance explained by the sum of squares
%
% Usage :
% =======
% [avgmat, anovamat] = sadecomp(X, Ypd, adecomp_model);
%
% Related functions :
% ===================
% adecomp.m (performs ANOVA decomposition)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
%
% =========================================================================

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : size(Ypd,2)
    Ypd(:,i) = yformatconv(Ypd(:,i),'intvec');
end

% Allocates space for the output arguments
anovamat = cell(1, length(adecomp_model.anovamat));
avgmat = cell(1, length(adecomp_model.anovamat)-1);

% Stores input matrix
anovamat{1} = X;

% Calculates sum of squares of the initial matrix
ssq{1} = sum(anovamat{1}(:).^2); % total sum of squares

% Number of decompositions (the first position in anovamat is the 
% undeflated original matrix and the last is the residuals matrix)
ndecomp = length(adecomp_model.anovamat)-1; 

% Loops over the number of decompositions to perform
for j = 1 : ndecomp
    
    if j == 1 % for grand mean
        
        avgmat{j} = repmat(adecomp_model.avgmat{1}(1,:),[size(Ypd,1),1]);
        anovamat{1,j+1} = X - avgmat{j};
        ssq{2} = sum(avgmat{1}(:).^2); % total sum of squares
        
    else % for all other factors/interactions
        
        avgmat{j} = zeros(size(X));
        for k = 1 : size(adecomp_model.Ysm{j-1}{1},1)
            idx = ismember(Ypd(:, unique(adecomp_model.matid{j})), adecomp_model.Ysm{j-1}{1}(k,:),'rows');
            avgmat{j}(idx,:) = repmat(adecomp_model.Ysm{j-1}{2}(k,:),[sum(idx==1),1]);
        end
        
        anovamat{1,j+1} = anovamat{1,j} - avgmat{j};
        ssq{1,j+1} = sum(avgmat{1,j}(:).^2);
        ssqvarexp{1,j-1} = (ssq{1,j+1} / (ssq{1,1}-ssq{1,2})) * 100;
        
    end
    
end

ssq{1,end+1} = sum(anovamat{1,end}(:).^2); 
ssqvarexp{1,end+1} = (ssq{1,end} / (ssq{1,1}-ssq{1,2})) * 100; 

end
