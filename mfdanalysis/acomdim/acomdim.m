function [acomdim_res] = acomdim(X, Y, Options)
%
% The acomdim function performs ANOVA-ComDim (AComdim) according to
% Jouan-Rimbaud Bouveresse et al. (2011). This method is based on APCA
% except for the step using PCA. In fact, AComDim treats all tables in
% one step instead of calculating a PCA model for each mean effect table to
% which the residuals matrix was added.
%
% AComDim analyzes in one step the main effect/interaction matrices along 
% with the residuals matrix.
% 
% References :
% ============
% Jouan-Rimbaud Bouveresse, D., Pinto, R. C., Schmidtke, L. M., 
% Locquet, N., & Rutledge, D. N. (2011). Identification of significant 
% factors by an extension of ANOVA-PCA based on multi-block analysis. 
% Chemometrics and Intelligent Laboratory Systems, 106(2), 173?182. 
% https://doi.org/10.1016/j.chemolab.2010.05.005
%
% de Figueiredo, M., Giannoukos, S., Wüthrich, C., Zenobi, R., & 
% Rutledge, D. N. (s. d.). A tutorial on the analysis of multifactorial 
% designs from one or more data sources using AComDim. Journal of 
% Chemometrics, n/a(n/a), e3384. https://doi.org/10.1002/cem.3384
%
% Input arguments :
% =================
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample 
%     identified by consecutive integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.ccs : the number of common components to extract
%   Options.decomp : ANOVA decomposition method 'classical', 'glm',
%   'rebalanced';
%   Options.coding : if Options.decomp is 'glm', the coding method of the
%   design matrices can be chosen as 'sumcod' or 'wecod';
%   Options.interactions : the maximum number of interactions to consider
%   Options.permtest : permutation tests for significance testing if ~= 0
%   Options.nperms : number of permutations if Options.permtest ~= 0
% 
% Output arguments :
% ==================
% acomdim_res : structure containing results of the AComdDim
%   acomdim_res.Xprepro : stores the preprocessed data before decomposition
%   acomdim_res.Y: numerical array of the experimental design;
%   acomdim_res.Y: numerical array of the experimental data;
%   acomdim_res.Xm: the grand mean matrix;
%   acomdim_res.Xf: cell array of the pure effect matrices;
%   acomdim_res.Xe: numerical array of the pure error (residuals);
%   acomdim_res.ssq: sum of squares of the effects;
%   acomdim_res.ssqvarexp: sum of squares explained variation;
%   acomdim_res.Xfaug: cell array of the augmented effect matrices;
%   acomdim_res.effects: cell array of the effects (factor indices);
%
% If the ANOVA decomposition uses the GLM methodology, also contains:
%   acomdim_res.Ef: the residuals without considering the effect f in the model;
%   acomdim_res.cod: the coding of the experimental design Y;
%   acomdim_res.B: the GLM parameters;
%
%   acomdim_res.scores : global scores of AComDim
%   acomdim_res.loadings : global loadings of AComDim
%   acomdim_res.lambda : saliences of AComDim
%   acomdim_res.maxlambda : explained effect by each salience
%   acomdim_res.varexp : explained variance of AComDim 
%   acomdim_res.F : F-ratios and p values for significance testing
%   acomdim_res.comdim : All ComDim results, see comdimc.m for more information
%   acomdim_res.Options : options used to perform AComDim
%
%   if Options.permtest == 1 :
%   acomdim_res.ptest : stores permutation tests results and p-values if
%       performed
%
%   acomdim_res.Options : options used to perform AComDim
% 
% Usage :
% =======
% AComDim : classical ANOVA decomposition
% Options.decomp = 'classical';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [acomdim_res] = acomdim(X, Y, Options);
%
% or 
%
% AComDim+ : sum-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [acomdim_res] = acomdim(X, Y, Options);
%
% or 
%
% WE-AComDim : weighted-effect-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'swecod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [acomdim_res] = acomdim(X, Y, Options);
%
% or
%
% RAComDim : Rebalanced AComDim ANOVA decomposition
% Options.decomp = 'rebalanced';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% [acomdim_res] = acomdim(X, Y, Options);
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

% Checks if number of ccs to extract was defined
if isfield(Options,'ccs') == 0
    Options.ccs = 9; % arbitrary
end

%% Preprocessing of the input data

[acomdim_res.Xprepro] = xpreproc(X, Options); % preprocessing as chosen

%% Performs ANOVA decomposition of each variable into effect matrices

% Decomposition is either classical, using GLM methodology or rebalancing procedure
if strcmp(Options.decomp,'rebalanced') == 1
    [adecomp_res] = radecomp(acomdim_res.Xprepro.data, Y, Options); % see radecomp.m for information
else
    [adecomp_res] = adecomp(acomdim_res.Xprepro.data, Y, Options); % see adecomp.m for information
end

% Copies fields of adecomp structure to the acomdim_res structure
fieldn = fieldnames(adecomp_res)';
for fn = fieldn(1:end-1) % without options
   acomdim_res.(fn{1}) = adecomp_res.(fn{1});
end

%% Performs ComDim on ANOVA matrices after adding residuals matrix 

neff = length(acomdim_res.effects);

% Extracts augmented effect matrices and adds matrix of residuals 
Xs = [acomdim_res.Xfaug,acomdim_res.Xe];

% Centers and normalizes tables to unit variance
for i = 1 : neff + 1 % +1 because residuals are taken into account
    cmdata = mean(Xs{i},1); % Column mean calculation
    Xs{i} = Xs{i} - ( ones(size(Xs{i},1),1) * cmdata );
    Xs{i} = Xs{i} ./ sqrt( sum( sum( Xs{i}.^2 ) ) ); % normalization of table i
end

% Performs ComDim on all tables simultaneously
[comdim_res] = comdimc(Xs, Options.ccs, Options);

acomdim_res.comdim = comdim_res; % ComDim full model
acomdim_res.scores = comdim_res.Q; % ComDim local scores
acomdim_res.loadings = comdim_res.P; % ComDim local loadings
acomdim_res.lambda = comdim_res.saliences; % ComDim block saliences
acomdim_res.varexp = comdim_res.varexp.Xsglobal; % CCs explained variance

%% Calculation of F from the saliences

maxsal = zeros(1,acomdim_res.comdim.CCs);
for i = 1 : acomdim_res.comdim.CCs
    maxsal(i) = find(acomdim_res.comdim.saliences(:,i) == max(acomdim_res.comdim.saliences(:,i)));
end

acomdim_res.maxlambda = maxsal; % Maximum saliences identifying the explained effects

% Uses all CCs linked to residuals
maxsalcc = find(maxsal == (neff+1));

Fnum = ones(neff+1,1) * sum(acomdim_res.comdim.saliences(neff+1,maxsalcc),2);
Fden = sum(acomdim_res.comdim.saliences(:,maxsalcc),2);

% Only uses CC1 to calculate the F-ratios
maxsalcc = 1; 

Fnumall = ones(neff+1,1) * sum(acomdim_res.comdim.saliences(neff+1,maxsalcc),2);
Fdenall = sum(acomdim_res.comdim.saliences(:,maxsalcc),2);

acomdim_res.F.val = Fnum ./ Fden;
acomdim_res.F.pval = 1 - fcdf(acomdim_res.F.val, n-1, n-1);
acomdim_res.F.valall = Fnumall ./ Fdenall;
acomdim_res.F.pvalall = 1 - fcdf(acomdim_res.F.valall, n-1, n-1);
acomdim_res.F.crit95 = finv(0.95, n-1, n-1);
acomdim_res.F.crit99 = finv(0.99, n-1, n-1);

acomdim_res.maxsal = maxsal;

% Note that the fields F.val/F.pval refer to the significance tests based
% on CC1 solely for the calculation of the F-ratios. 
% F.valall/F.pvalall refer to the use of all CCs linked to the residuals.

%% In case Rebalanced AComDim was used, the algorithm enters this section

% If the rebalanced method was used, the augmented effect matrices are
% projected on the ComDim model
if strcmp(Options.decomp,'rebalanced') == 1
    
    % ANOVA decomposition of the original matrix according to the mean
    % effects of each factor and interaction
    [sadecomp_res] = sadecomp(acomdim_res.Xprepro.data, Y, adecomp_res);

    % Copies fields of adecomp structure to the asca_res structure
    fieldn = fieldnames(sadecomp_res)';
    for fn = fieldn(1:end-1) % without options
        acomdim_res.(fn{1}) = sadecomp_res.(fn{1});
    end

    % Augmented effect matrices + matrix of residuals
    Xs = [sadecomp_res.Xfaug,sadecomp_res.Xe]; % projection

    % Projects initial data into ComDim model space
    [comdim_proj] = comdimp(Xs, comdim_res);
    acomdim_res.scores = comdim_proj.Q; % updates the scores with those of the initial matrix

end

acomdim_res.Options = Options;

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    
    [acomdim_res.ptest] = mfdptest(acomdim_res);
    
end

end