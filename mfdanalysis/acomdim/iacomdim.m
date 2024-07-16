function [iacomdim_res] = iacomdim(Xs, Yrep, Ys, Options)
%
% The iacomdim function performs ANOVA-ComDim (AComdim) according to
% de Figueiredo et al. (2022). This method is based on APCA
% except for the step using PCA. In fact, AComDim treats all tables in
% one step instead of calculating a PCA model for each mean effect table to
% which the residuals matrix was added.
%
% AComDim analyzes in one step the main effect/interaction matrices along 
% with the residuals matrix.
%
% What differentiates the iacomdim function from the acomdim is that it
% is an extension to cases where data from multiple sources are available.
% 
% References :
% ============
% de Figueiredo, M., Giannoukos, S., Wüthrich, C., Zenobi, R., & 
% Rutledge, D. N. (s. d.). A tutorial on the analysis of multifactorial 
% designs from one or more data sources using AComDim. Journal of 
% Chemometrics, n/a(n/a), e3384. https://doi.org/10.1002/cem.3384
%
% Input arguments :
% =================
% Xs : cell array of data matrices with samples in the rows and variables in the columns
%     where each cell contains one source of data
%
% Ys : cell array of the design of experiments associated with each table
% in Xs
%
% Yrep : cell array of integer vectors identifying replicates accross all
% tables
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
% Options.prepro.X.type = {{'nothing'}}; % or 'cmeancenter' or 'autoscale' 
% Options.decomp = 'classical'; % or 'glm'
% Options.interactions = 2; % number of interactions to consider
% Options.permtest = 1; % performs permutations tests (1) or not (0)
% Options.nperms = 10000; % number of permutations if Options.permtest == 1
% Options.ccs = 10; % number of common components to extract
% [iacomdim_res] = iacomdim(Xs, Y, Options);
%
% Related function :
% ==================
% xpreproc.m (performs preprocessing)
% acomdim.m (performs AComDim on a single source)
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

q = size(Ys{1},2);

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if Yrep was provided, otherwise assumes aligned replicates in thetables
if isempty(Yrep)
    Yrep = repmat({(1:size(Xs{1}))'},[1,length(Xs)]);
end

% Check is Ys is provided as a numeric array, if so assumes the same design
% for all tables
if isnumeric(Ys)
    Ys = repmat({Ys},[1,length(Xs)]);
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

% Checks if number of ccs to extract was defined
if isfield(Options,'ccs') == 0
    Options.ccs = 9; % arbitrary
end

%% Preprocessing of the input data

[iacomdim_res.Xprepro] = xpreproc(Xs, Options); % preprocessing as chosen

%% Performs multivariate ANOVA decomposition into effect matrices

% Decomposition is either classical or uses the GLM methodology

if strcmp(Options.decomp,'rebalanced') == 1
    [adecomp_res, Xs, Ys] = iradecomp({iacomdim_res.Xprepro.data}, Yrep, Ys, Options); % see iradecomp.m for information
    Y = Ys{1};
else
    for i = 1 : length(Xs)
        [adecomp_res(i)] = adecomp(iacomdim_res.Xprepro(i).data, Ys{i}, Options); % see adecomp.m for information
    end
    Y = Ys{1};
end

iacomdim_res.X = Xs;
iacomdim_res.Y = Y;
iacomdim_res.adecomp = adecomp_res;

%% Performs ComDim on anova matrices after adding residuals matrix + the residuals block

% Concatenates all augmented effect matrices and residuals of all tables
Xsall = {adecomp_res.Xfaug}; 
for i = 1 : length(Xsall); Xsall{i} = [Xsall{i},adecomp_res(i).Xe]; end
Xsall = cat(2,Xsall{:});

for i = 1 : length(Xsall)
    cmdata = mean(Xsall{i},1); % Column mean calculation
    Xsall{i} = Xsall{i} - ( ones(size(Xsall{i},1),1) * cmdata );
    Xsall{i} = Xsall{i} ./ sqrt( sum( sum( Xsall{i}.^2 ) ) ); % normalization of table i
end

% Performs ComDim
[comdim_res] = comdimc(Xsall, Options.ccs, Options);

iacomdim_res.comdim = comdim_res; % ComDim full model
iacomdim_res.scores = comdim_res.Q; % ComDim local scores
iacomdim_res.loadings = comdim_res.P; % ComDim local loadings
iacomdim_res.lambda = comdim_res.saliences; % ComDim block saliences
iacomdim_res.varexp = comdim_res.varexp.Xsglobal; % CCs explained variance

%% Calculation of F from the saliences

matid = [iacomdim_res.adecomp(1).effects,'Residuals'];
matid2 = repmat(1:length(matid),[1,size(iacomdim_res.adecomp,2)]);

maxsal = zeros(1,iacomdim_res.comdim.CCs);
for i = 1 : iacomdim_res.comdim.CCs
    idx = find(iacomdim_res.comdim.saliences(:,i) == max(iacomdim_res.comdim.saliences(:,i)));
    maxsal(i) = idx(1);
    maxsal(i) = matid2(maxsal(i));
end

iacomdim_res.maxlambda = maxsal; % Maximum saliences identifying the explained effects

% Uses all CCs linked to residuals
maxsalcc = find(maxsal == length(matid));

Fnum = reshape((sum(iacomdim_res.comdim.saliences(matid2 == length(matid),maxsalcc),2) * ones(length(matid),1)')',[length(matid)*size(iacomdim_res.adecomp,2),1]);
Fden = sum(iacomdim_res.comdim.saliences(:,maxsalcc),2);

% Only uses CC1 to calculate the F-ratios
maxsalcc = 1;

Fnumall = reshape((sum(iacomdim_res.comdim.saliences(matid2 == length(matid),maxsalcc),2) * ones(length(matid),1)')',[length(matid)*size(iacomdim_res.adecomp,2),1]);
Fdenall = sum(iacomdim_res.comdim.saliences(:,maxsalcc),2);

n = size(iacomdim_res.adecomp(1).X,1);
iacomdim_res.F.id = [cellstr(string(matid2));matid(matid2)]';
iacomdim_res.F.val = Fnum ./ Fden;
iacomdim_res.F.valall = Fnumall ./ Fdenall;
iacomdim_res.F.pval = 1 - fcdf(iacomdim_res.F.val, n-1, n-1);
iacomdim_res.F.pvalall = 1 - fcdf(iacomdim_res.F.valall, n-1, n-1);
iacomdim_res.F.crit95 = finv(0.95, size(Y,1)-1, size(Y,1)-1);
iacomdim_res.F.crit99 = finv(0.99, size(Y,1)-1, size(Y,1)-1);

% Note that the fields F.val/F.pval refer to the significance tests based
% on CC1 solely for the calculation of the F-ratios. 
% F.valall/F.pvalall refer to the use of all CCs linked to the residuals.

%% In case Rebalanced AComDim was used, the algorithm enters this section

% If the rebalanced method was used, the augmented effect matrices are
% projected on the ComDim model
if strcmp(Options.decomp,'rebalanced') == 1
    
    % ANOVA decomposition of the original matrix according to the mean
    % effects of each factor and interaction
    [sadecomp_res] = isadecomp(iacomdim_res.X, Y, adecomp_res);
    iacomdim_res.adecomp = sadecomp_res;
    [iacomdim_res.adecomp.nans] = adecomp_res.nans;

    % Augmented effect matrices + matrix of residuals
    Xsall = {sadecomp_res.Xfaug}; 
    for i = 1 : length(Xsall); Xsall{i} = [Xsall{i},sadecomp_res(i).Xe]; end
    Xsall = cat(2,Xsall{:});

    % Projects initial data into ComDim model space
    [comdim_proj] = comdimp(Xsall, iacomdim_res.comdim);
    iacomdim_res.scores = comdim_proj.Q; % updates the scores with those of the initial matrix

end

iacomdim_res.Options = Options;

%% Reshapes the saliences and the loadings for more clarity

ntabs = length(adecomp_res(1).effects)+1; % effects + residuals
nblocks = length(Xs); % number of sources

% Lambda weights (saliences)
iacomdim_res.lambda = permute(reshape(iacomdim_res.lambda,[ntabs,nblocks,Options.ccs]),[2,1,3]);

% We rely on the fact that the global loadings are the concatenated local ones
iacomdim_res.loadings = reshape(iacomdim_res.comdim.Ploc,[ntabs,nblocks])';

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    
    for i = 1 : length(Xs)
        tmp = iacomdim_res.adecomp(i); 
        tmp.Y = iacomdim_res.Y;
        tmp.Options = iacomdim_res.Options;
        [iacomdim_res.ptest(i)] = mfdptest(tmp);
    end
    
end


end