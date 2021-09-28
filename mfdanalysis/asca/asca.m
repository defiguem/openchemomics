function [asca_res] = asca(X, Y, Options)
%
% The asca function performs ANOVA-Simultaneous Components Analysis (ASCA).
% In this implementation, an ANOVA decomposition into main factors and
% interactions is performed. Then, in the SCA part, a Singular Values 
% Decomposition (SVD) is applied to the effect matrices with and without
% prior addition of residuals.
%
% The authors first implemented ASCA using the scores directly calculated
% from the SVD on the effect matrices. A second way to calculate these 
% scores is by first adding the residuals matrix to the effect matrices,
% generating augmented effect matrices, then using the previously 
% calculated loadings to calculate the new scores. These scores restitute
% the natural variability of the data.
%
% The ANOVA decomposition step can be peformed using different strategies.
% For more information on the methods available see related functions
% below.
% 
% References :
% ============
% Smilde, A. K., Jansen, J. J., Hoefsloot, H. C. J., Lamers, R.-J. A. N., 
% van der Greef, J., & Timmerman, M. E. (2005). ANOVA-simultaneous 
% component analysis (ASCA): A new tool for analyzing designed 
% metabolomics data. Bioinformatics, 21(13), 3043?3048. 
% https://doi.org/10.1093/bioinformatics/bti476
% 
% Jansen, J. J., Hoefsloot, H. C. J., van der Greef, J., Timmerman, M. E., 
% Westerhuis, J. A., & Smilde, A. K. (2005). ASCA : Analysis of multivariate 
% data obtained from an experimental design. Journal of Chemometrics, 19(9), 
% 469?481. https://doi.org/10.1002/cem.952
%
% Vis, D. J., Westerhuis, J. A., Smilde, A. K., & van der Greef, J. (2007). 
% Statistical validation of megavariate effects in ASCA. BMC Bioinformatics, 
% 8(1). https://doi.org/10.1186/1471-2105-8-322
%
% Thiel, M., Féraud, B., & Govaerts, B. (2017). ASCA+ and APCA+: 
% Extensions of ASCA and APCA in the analysis of unbalanced multifactorial 
% designs: Analyzing unbalanced multifactorial designs with ASCA+ and APCA+. 
% Journal of Chemometrics, 31(6), e2895. https://doi.org/10.1002/cem.2895
%
% Ali, N., Jansen, J., van den Doel, A., Tinnevelt, G. H., & Bocklitz, T. 
% (2020). WE-ASCA?: The Weighted-Effect ASCA for Analyzing Unbalanced 
% Multifactorial Designs?A Raman Spectra-Based Example. Molecules, 26(1). 
% https://doi.org/10.3390/molecules26010066
%
% Input arguments :
% =================
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample 
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.decomp : ANOVA decomposition method 'classical', 'glm',
%   'res_classical', 'res_glm';
%   Options.coding : if Options.decomp is 'glm', the coding method of the
%   design matrices can be chosen as 'sumcod' or 'wecod';
%   Options.interactions : the maximum number of interactions to consider
%   Options.permtest : permutation tests for effects significance if ~= 0;
%   Options.nperms : number of permutations to perform 
%   if Options.permtest == 1;
%   Options.nres : number of resamplings to perform if Options.decomp is
%   'res_classical' or 'res_glm';
% 
% Output arguments :
% ==================
% asca_res : structure containing results of the ASCA 
%   asca_res.Y : stores the design matrix
%   asca_res.Xprepro : stores the preprocessed data before decomposition
%   asca_res.anovamat : ANOVA matrices contained in cells (the last cell
%       contains the residuals)
%   asca_res.avgmat : cell structure of factor/interation effects
%       (see asca_res.matid for information on factor/interactions)
%   asca_res.matid : identity of ANOVA matrices (grand mean, effects,
%       interactions and residual pure error)
%   asca_res.ssq : sum of squares for each effect based on asca_res.avgmat
%   asca_res.ssqvarexp : variance explained by the sum of squares
%   asca_res.Ys : stores the design matrices of factors and interactions
%   asca_res.Xs : stores the effect matrix to which the pure error 
%       was added
%   asca_res.Ysm : mean effect vector for each factor and interaction
%   asca_res.scores : scores of SCA for ANOVA matrices without adding  
%       pure error and performing SCA
%   asca_res.resscores : scores of SCA for ANOVA matrices after adding 
%       pure error and performing SCA
%   asca_res.loadings : loadings of SCA for ANOVA matrices without adding  
%       pure error and performing SCA
%   asca_res.evals : eigenvalues of SCA for ANOVA matrices without adding  
%       pure error and performing SCA
%   asca_res.varexp : explained variance of SCA for ANOVA matrices without 
%       adding   pure error and performing SCA
%   asca_res.Options : options used to perform ASCA
%
%   if Options.permtest == 1 :
%   asca_res.ptest : stores permutation tests results and p-values if
%       performed
%
%   if Options.decomp == 'res_classical' || 'res_glm':
%   asca_res.pavgmat : cell structure of factor/interation effects
%   asca_res.panovamat : deflated ANOVA matrices contained in cells
%   asca_res.pssq : sum of squares for each effect based on the synthetic
%   effect matrices
%   asca_res.pssqvarexp : variance explained by the sum of squares of the
%   synthetic effect matrices
%   asca_res.pscores : scores of mean resampled data matrix 
%   asca_res.presscores : scores of original data matrix after adding back
%   the matrix of residuals
% 
% Usage :
% =======
% ASCA : classical ANOVA decomposition
% Options.decomp = 'classical';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [asca_res] = asca(X, Y, Options);
%
% or 
%
% ASCA+ : sum-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [asca_res] = asca(X, Y, Options);
%
% or 
%
% WE-ASCA : weighted-effect-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [asca_res] = asca(X, Y, Options);
%
% or
%
% RASCA : Resampled ASCA classical ANOVA decomposition
% Options.decomp = 'res_classical';
% Options.nres = 100; % number of resamplings to perform
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [asca_res] = asca(X, Y, Options);
%
% Related function :
% ==================
% xpreproc.m (preprocessing of the input data matrix)
% adecomp.m (performs ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
% radecomp.m (performs multiple resamplings before ANOVA decomposition)
% sadecomp.m (ANOVA decomposition of the orignal matrix for Resampled ASCA)
% mfdptest.m (performs permutation tests)
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

[n,p] = size(X); % size of X
[~,q] = size(Y); % size of Y

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if the multivariate anova decomposition method was defined
if isfield(Options,'decomp') == 0
    Options.decomp = 'classical';
end

% Checks if the maximum number of interactions to consider was defined
if isfield(Options,'interactions') == 0
    Options.interactions = 2;
elseif isfield(Options,'interactions') == 1 && Options.interactions > q
    Options.interactions = q; % maximum number of interactions cannot be > q
end

% Checks if permutation tests must be conducted
if isfield(Options,'permtest') == 0
    Options.permtest = 0;
end

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : q
    [res] = classcol(Y(:,i));
    Y(:,i) = res.classvec;
end

asca_res.Y = Y; % stores the design of the experiments

%% Preprocessing of the input data

[asca_res.Xprepro] = xpreproc(X, Options); % preprocessing as chosen

%% Performs multivariate ANOVA decomposition into effect matrices

% Decomposition is either classical, using GLM methodology or using a
% resampling procedure
if strcmp(Options.decomp,'res_classical') == 1 || strcmp(Options.decomp,'res_glm') == 1
    [adecomp_res] = radecomp(asca_res.Xprepro.data, Y, Options); % see radecomp.m for information
else
    [adecomp_res] = adecomp(asca_res.Xprepro.data, Y, Options); % see adecomp.m for information
end

% Copies fields of adecomp structure to the asca_res structure
for fn = fieldnames(adecomp_res)'
   asca_res.(fn{1}) = adecomp_res.(fn{1});
end

%% Performs SCA by SVD on each ANOVA matrix after adding residuals matrix

ntab = length(asca_res.avgmat) - 1; % note that SCA is not done on the first table corresponding to the grand mean matrix

for i = 1 : ntab

    % Performs SCA
    R = rank(asca_res.avgmat{1,i+1}); % rank of the average matrix
    [U,S,P] = svds(asca_res.avgmat{1,i+1}, R); % svd decomposition
    T = U * S; % scores
    resT = (asca_res.avgmat{1,i+1} + asca_res.anovamat{1,end}) * P; % scores after adding back the residuals' matrix
    evals = diag( S ) ./ (n-1); % eigenvalues
    varexp = 100 * ( evals ./ sum(evals) ); % explained variance
    
    % Stores residuals augmented scores, scores, loadings and explained variance
    asca_res.scores{1,i} = T; % scores without residuals
    asca_res.resscores{1,i} = resT; % scores with residuals
    asca_res.loadings{1,i} = P; % loadings without residuals
    asca_res.evals{1,i} = evals; % eigenvalues
    asca_res.varexp{1,i} = varexp; % explained variance
    
end

asca_res.Options = Options; % stores options used

% If the resampling method was used, the scores of the initial matrix are
% projected on all SCA models
if strcmp(Options.decomp,'res_classical') == 1 || strcmp(Options.decomp,'res_glm') == 1
    
    % ANOVA decomposition of the original matrix according to the mean
    % effects of each factor and interaction
    [asca_res.pavgmat, asca_res.panovamat, asca_res.pssq, asca_res.pssqvarexp] = sadecomp(asca_res.Xprepro.data, Y, adecomp_res, Options);
    
    ntab = length(asca_res.avgmat)-1; % note that pca is not done on the firts table corresponding to the grand mean matrix
    
    for i = 1 : ntab
        
        % Projects initial matrix in the SCA model
        resTp = asca_res.pavgmat{1,i+1} * asca_res.loadings{1,i}; % note that SCA is not done on the first table corresponding to the grand mean matrix
        
        % Stores residuals augmented scores of original matrix
        asca_res.pscores{1,i} = resTp;
        
        % Projects initial matrix in the SCA model
        resTp = (asca_res.pavgmat{1,i+1} + asca_res.panovamat{1,end}) * asca_res.loadings{1,i}; % note that SCA is not done on the first table corresponding to the grand mean matrix
        
        % Stores residuals augmented scores of original matrix
        asca_res.presscores{1,i} = resTp;
        
    end
end

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    [asca_res.ptest] = mfdptest(asca_res);
end

end