function [apca_res] = apca(X, Y, Options)
%
% The apca function performs ANOVA-Principal Components Analysis (APCA).
% From the authors (Harrington et al., 2005) : 
% "Covariations are separated using ANOVA into main effects and 
% interaction. The covariances for each effect are combined with the pure 
% error and subjected to PCA. If the main effect is significant compared to 
% the residual error, the first principal component will span this source 
% of variation".
%
% PCA is performed using Singular Values Decomposition (SVD).
%
% Note that in this implementation, the General Linear Model (GLM) can be 
% used for decomposition instead of the classical ANOVA decomposition. 
% This approach is based on the paper of Thiel et al. (2017) and has the 
% advantage of providing unbiased estimators even with unbalanced designs. 
% In the cases where the design is balanced, these estimators are identical 
% to those given by the classical ANOVA decomposition. Read reference for 
% more information.
% 
% Reference :
% ===========
% Harrington, P. de B., Vieira, N. E., Espinoza, J., Nien, J. K., Romero, 
% R., & Yergey, A. L. (2005). Analysis of variance?principal component 
% analysis: A soft tool for proteomic discovery. Analytica Chimica Acta, 
% 544(1-2), 118-127. https://doi.org/10.1016/j.aca.2005.02.042
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
% de Figueiredo, M., Giannoukos, S., Rudaz, S., Zenobi, R., & Boccard, J. 
% (s. d.). Efficiently handling high-dimensional data from multifactorial 
% designs with unequal group sizes using Rebalanced ASCA (RASCA). 
% Journal of Chemometrics, n/a(n/a), e3401. https://doi.org/10.1002/cem.3401
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
%   'rebalanced';
%   Options.coding : if Options.decomp is 'glm', the coding method of the
%   design matrices can be chosen as 'sumcod' or 'wecod';
%   Options.interactions : the maximum number of interactions to consider
%   Options.permtest : permutation tests for effects significance if ~= 0;
%   Options.nperms : number of permutations to perform 
%
% Output arguments :
% ==================
% apca_res : structure containing results of the APCA
%   apca_res.Xprepro : stores the preprocessed data before decomposition
%   apca_res.Y: numerical array of the experimental design;
%   apca_res.Y: numerical array of the experimental data;
%   apca_res.Xm: the grand mean matrix;
%   apca_res.Xf: cell array of the pure effect matrices;
%   apca_res.Xe: numerical array of the pure error (residuals);
%   apca_res.ssq: sum of squares of the effects;
%   apca_res.ssqvarexp: sum of squares explained variation;
%   apca_res.Xfaug: cell array of the augmented effect matrices;
%   apca_res.effects: cell array of the effects (factor indices);
%
% If the ANOVA decomposition uses the GLM methodology, also contains:
%   apca_res.Ef: the residuals without considering the effect f in the model;
%   apca_res.cod: the coding of the experimental design Y;
%   apca_res.B: the GLM parameters;
%
%   apca_res.scores : scores of PCA for ANOVA matrices 
%   apca_res.loadings : loadings of PCA for ANOVA matrices 
%   apca_res.evals : eigenvalues of PCA for ANOVA matrices 
%   apca_res.varexp : explained variance of PCA for ANOVA matrices 
%   apca_res.Options : options used to perform APCA
%
%   if Options.permtest == 1 :
%   apca_res.ptest : stores permutation tests results and p-values if
%       performed
% 
% Usage :
% =======
% APCA : classical ANOVA decomposition
% Options.decomp = 'classical';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [apca_res] = apca(X, Y, Options);
%
% or 
%
% APCA+ : sum-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [apca_res] = apca(X, Y, Options);
%
% or 
%
% WE-APCA : weighted-effect-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'swecod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [apca_res] = apca(X, Y, Options);
%
% or
%
% RAPCA : Rebalanced APCA ANOVA decomposition
% Options.decomp = 'rebalanced';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [apca_res] = apca(X, Y, Options);
%
% Related function :
% ==================
% xpreproc.m (preprocessing of the input data matrix)
% adecomp.m (performs ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
% radecomp.m (performs rebalanced ANOVA decomposition)
% sadecomp.m (performs synthetic ANOVA decomposition)
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
    Options.decomp = 'rebalanced';
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

%% Preprocessing of the input data

[apca_res.Xprepro] = xpreproc(X, Options); % preprocessing as chosen

%% Performs ANOVA decomposition of each variable into effect matrices

% Decomposition is either classical, using GLM methodology or rebalancing procedure
if strcmp(Options.decomp,'rebalanced') == 1
    [adecomp_res] = radecomp(apca_res.Xprepro.data, Y, Options); % see radecomp.m for information
else
    [adecomp_res] = adecomp(apca_res.Xprepro.data, Y, Options); % see adecomp.m for information
end

% Copies fields of adecomp structure to the apca_res structure
fieldn = fieldnames(adecomp_res)';
for fn = fieldn(1:end-1) % without options
   apca_res.(fn{1}) = adecomp_res.(fn{1});
end

%% Performs PCA by SVD on each anova matrix after adding residuals matrix

neff = length(apca_res.effects); % note that pca is not done on the firts table corresponding to the grand mean matrix
PCs = min([n,p]); % automated to min([n,p]

for i = 1 : neff
    
    % Performs PCA by SVD
    Options.pcatype = 'svd'; 
    [pca_model] = pcac(apca_res.Xf{i} + apca_res.Xe, PCs, Options); % adds residuals to effect matrices before PCA
    
    % Stores scores, loadings and explained variance
    apca_res.scores{1,i} = pca_model.scores;
    apca_res.loadings{1,i} = pca_model.loadings;
    apca_res.evals{1,i} = pca_model.evals;
    apca_res.varexp{1,i} = pca_model.varexp;
    
end

apca_res.Options = Options; % stores options used

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    [apca_res.ptest] = mfdptest(apca_res);
end

end