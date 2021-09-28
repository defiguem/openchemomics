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

% Switches cases depending on the multivariate ANOVA decomposition method
switch Options.decomp
    
    case {'classical','res_classical'}
        
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
        
    case {'glm','res_glm'}
        
        % Either codes design matrices using the sum coding or the
        % weighted-effect coding
        if strcmp(Options.coding,'sumcod') == 1
            
            % Sum coding of the input design matrix according to 
            % Thiel et al. (2017)
            [cod, codempty] = sumcoding(Ypd, Options);
            
        elseif strcmp(Options.coding,'wecod') == 1
            
            % Weighted-effect coding of the input design matrix according
            % to Nieuwenhuis et al. (2017)
            [cod, codempty] = wecoding(Ypd, Options);
            
        end
        
        % Creates empty B (it is just useful for later calculations...)
        bempty = {};
        for i = 1 : length(cod)
            bempty{i} = zeros(size(X,2),size(cod{i},2));
        end
        
        % Loops over the number of decompositions to perform
        for j = 1 : ndecomp
            
            if j == 1 % for grand mean
                
                % Grand mean matrix : initialization of decomposition into GLM matrices
                B = bempty;
                B{1} = adecomp_model.B{1};
                Xf = codempty;
                Xf{1} = cod{1};
                
                avgmat{j} = cat(2,Xf{:}) * cat(2,B{:})'; % grand mean matrix 
                anovamat{1,j+1} = X - avgmat{j};
                ssq{2} = sum(avgmat{1}(:).^2); 
                
            else % for all other factors/interactions
                
                B = bempty;
                B{j} = adecomp_model.B{j};
                Xf = codempty;
                Xf{j} = cod{j};
            
                avgmat{j} = cat(2,Xf{:}) * cat(2,B{:})';
                anovamat{j+1} = anovamat{1,j} - avgmat{j};
                ssq{1,j+1} = sum(avgmat{1,j}(:).^2); 
                ssqvarexp{1,j-1} = (ssq{1,j+1} / (ssq{1,1}-ssq{1,2})) * 100; 
                
            end
            
        end
        
end

ssq{1,end+1} = sum(anovamat{1,end}(:).^2); 
ssqvarexp{1,end+1} = (ssq{1,end} / (ssq{1,1}-ssq{1,2})) * 100; 

end
