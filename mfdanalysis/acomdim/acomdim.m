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
%   Options.decomp : 'classical' or 'glm', see function adecomp.m
%   Options.interactions : the maximum number of interactions to consider
%   Options.permtest : permutation tests for significance testing if ~= 0
%   Options.nperms : number of permutations if Options.permtest ~= 0
% 
% Output arguments :
% ==================
% acomdim_res : structure containing results of the A-ComDim 
%   acomdim_res.Y : stores the design matrix
%   acomdim_res.Xprepro : stores the preprocessed data before decomposition
%   acomdim_res.anovamat : progressively deflated ANOVA matrices contained 
%       in cells (the last cell contains the residuals)
%   acomdim_res.avgmat : cell structure of mean factor/interation effects
%       (see acomdim_res.matid for information on factor/interactions)
%   acomdim_res.matid : identity of ANOVA matrices (grand mean, effects,
%       interactions and residual pure error)
%   acomdim_res.ssq : sum of squares for each effect based on acomdim_res.avgmat
%   acomdim_res.ssqvarexp : variance explained by the sum of squares
%   acomdim_res.Ys : stores the design matrices of factors and interactions
%   acomdim_res.Xs : stores the effect matrix to which the pure error 
%       was added
%   acomdim_res.Ysm : mean effect vector for each factor and interaction
%   acomdim_res.comdim : ComDim results (saliences, scores, loadings,
%   etc.), see comdimc.m for more information
%   acomdim_res.F : F-ratios and p values for significance testing
%   acomdim_res.Options : options used to perform AComDim
%   acomdim_res.ptest : stores permutation test results and p-values if
%       performed
% 
% Usage :
% =======
% Options.prepro.X.type = {{'nothing'}}; % or 'cmeancenter' or 'autoscale' 
% Options.decomp = 'classical'; % or 'glm'
% Options.interactions = 2; % number of interactions to consider
% Options.permtest = 1; % performs permutations tests (1) or not (0)
% Options.nperms = 10000; % number of permutations if Options.permtest == 1
% Options.ccs = 10; % number of common components to extract
% [acomdim_res] = acomdim(X, Y, Options);
%
% Related function :
% ==================
% xpreproc.m (performs preprocessing)
% adecomp.m (performs multivariate ANOVA decompostion)
% comdimc.m (performs ComDim)
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

% Checks if number of ccs to extract was defined
if isfield(Options,'ccs') == 0
    Options.ccs = 9; % arbitrary
end

% Checks if permutation tests must be conducted
if isfield(Options,'permtest') == 0
    Options.permtest = 0;
end

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : q
    Y(:,i) = yformatconv(Y(:,i),'intvec');
end

acomdim_res.Y = Y; % stores the design of the experiments

%% Preprocessing of the input data

[acomdim_res.Xprepro] = xpreproc(X, Options); % preprocessing as chosen

%% Performs multivariate ANOVA decomposition into effect matrices

% ANOVA decomposition (see adecomp.m for information)
[adecomp_res] = adecomp(acomdim_res.Xprepro.data, Y, Options);

% Copies fields of adecomp structure to the acomdim_res structure
for fn = fieldnames(adecomp_res)'
   acomdim_res.(fn{1}) = adecomp_res.(fn{1});
end

%% Performs ComDim on ANOVA matrices after adding residuals matrix 

ntab = length(acomdim_res.avgmat) - 1; % -1 because does not use grand mean matrix

% Extracts effects/interactions where last table is the residuals
Xs = adecomp_res.Xs;

% Centers and normalizes tables to unit variance
for i = 1 : ntab + 1 % +1 because residuals are taken into account
    cmdata = mean(Xs{i},1); % Column mean calculation
    Xs{i} = Xs{i} - ( ones(size(Xs{i},1),1) * cmdata );
    Xs{i} = Xs{i} ./ sqrt( sum( sum( Xs{i}.^2 ) ) ); % normalization of table i
end

% Performs ComDim on all tables simultaneously
[acomdim_res.comdim] = comdimc(Xs, Options.ccs, Options);

%% Calculation of F from the saliences

% The following 6 lines consider all CCs for which the block of residuals
% has the highest salience (otherwise it uses only CC1)
maxsal = zeros(1,acomdim_res.comdim.CCs);
for i = 1 : acomdim_res.comdim.CCs
    maxsal(i) = find(acomdim_res.comdim.saliences(:,i) == max(acomdim_res.comdim.saliences(:,i)));
end

% Uses all CCs linked to residuals
maxsalcc = find(maxsal == length(acomdim_res.matid)-1);

Fnum = ones(ntab+1,1) * sum(acomdim_res.comdim.saliences(ntab+1,maxsalcc),2);
Fden = sum(acomdim_res.comdim.saliences(:,maxsalcc),2);

% Only uses CC1 to calculate the F-ratios
maxsalcc = 1; 

Fnumall = ones(ntab+1,1) * sum(acomdim_res.comdim.saliences(ntab+1,maxsalcc),2);
Fdenall = sum(acomdim_res.comdim.saliences(:,maxsalcc),2);

acomdim_res.F.val = Fnum ./ Fden;
acomdim_res.F.pval = 1 - fcdf(acomdim_res.F.val, n-1, n-1);
acomdim_res.F.valall = Fnumall ./ Fdenall;
acomdim_res.F.pvalall = 1 - fcdf(acomdim_res.F.valall, n-1, n-1);
acomdim_res.F.crit95 = finv(0.95, n-1, n-1);
acomdim_res.F.crit99 = finv(0.99, n-1, n-1);

acomdim_res.Options = Options;

% Note that the fields F.val/F.pval refer to the significance tests based
% on CC1 solely for the calculation of the F-ratios. 
% F.valall/F.pvalall refer to the use of all CCs linked to the residuals.

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    
    [acomdim_res.ptest] = mfdptest(acomdim_res);
    
end

end