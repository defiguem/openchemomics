function [mbacomdim_res] = mbacomdim(Xs, Y, Options)
%
% The mbacomdim function performs ANOVA-ComDim (AComdim) according to
% Jouan-Rimbaud Bouveresse et al. (2011). This method is based on APCA
% except for the step using PCA. In fact, AComDim treats all tables in
% one step instead of calculating a PCA model for each mean effect table to
% which the residuals matrix was added.
%
% AComDim analyzes in one step the main effect/interaction matrices along 
% with the residuals matrix.
%
% What differentiates the mbacomdim function from the acomdim is that it
% is an extension to cases where data from multiple sources are available.
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
% Xs : cell array of data matrices with samples in the rows and variables in the columns
%     where each cell contains one source of data
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
%   acomdim_res.adecomp : structure with the ANOVA decomposition results of
%   each source of data, see adecomp.m for more information
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
% [mbacomdim_res] = mbacomdim(Xs, Y, Options);
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

[n,~] = size(Xs{1}); % size of X
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

mbacomdim_res.Y = Y; % stores the design of the experiments

%% Preprocessing of the input data

[mbacomdim_res.Xprepro] = xpreproc(Xs, Options); % preprocessing as chosen

%% Performs multivariate ANOVA decomposition into effect matrices

% Decomposition is either classical or uses the GLM methodology
for i = 1 : length(Xs)
    [adecomp_res(i)] = adecomp(mbacomdim_res.Xprepro(i).data, Y, Options); % see adecomp.m for information
end

mbacomdim_res.adecomp = adecomp_res;

%% Performs ComDim on anova matrices after adding residuals matrix + the residuals block

ntab = length(mbacomdim_res.adecomp(1).matid) - 1; % -1 because does not use grand mean matrix

Xsall = {adecomp_res.Xs}; 
Xsall = cat(2,Xsall{:});

for i = 1 : length(Xsall)
    cmdata = mean(Xsall{i},1); % Column mean calculation
    Xsall{i} = Xsall{i} - ( ones(size(Xsall{i},1),1) * cmdata );
    Xsall{i} = Xsall{i} ./ sqrt( sum( sum( Xsall{i}.^2 ) ) ); % normalization of table i
end

% Performs ComDim
[mbacomdim_res.comdim] = comdimc(Xsall, Options.ccs, Options);

%% Calculation of F from the saliences

matid = mbacomdim_res.adecomp(1).matid(2:end);
matid2 = repmat(1:length(matid),[1,size(mbacomdim_res.adecomp,2)]);

maxsal = zeros(1,mbacomdim_res.comdim.CCs);
for i = 1 : mbacomdim_res.comdim.CCs
    maxsal(i) = find(mbacomdim_res.comdim.saliences(:,i) == max(mbacomdim_res.comdim.saliences(:,i)));
    maxsal(i) = matid2(maxsal(i));
end

% Uses all CCs linked to residuals
maxsalcc = find(maxsal == length(matid));

Fnum = reshape((sum(mbacomdim_res.comdim.saliences(matid2 == length(matid),maxsalcc),2) * ones(length(matid),1)')',[length(matid)*size(mbacomdim_res.adecomp,2),1]);
Fden = sum(mbacomdim_res.comdim.saliences(:,maxsalcc),2);

% Only uses CC1 to calculate the F-ratios
maxsalcc = 1;

Fnumall = reshape((sum(mbacomdim_res.comdim.saliences(matid2 == length(matid),maxsalcc),2) * ones(length(matid),1)')',[length(matid)*size(mbacomdim_res.adecomp,2),1]);
Fdenall = sum(mbacomdim_res.comdim.saliences(:,maxsalcc),2);

mbacomdim_res.F.id = [cellstr(string(matid2));matid(matid2)]';
mbacomdim_res.F.val = Fnum ./ Fden;
mbacomdim_res.F.valall = Fnumall ./ Fdenall;
mbacomdim_res.F.pval = 1 - fcdf(mbacomdim_res.F.val, n-1, n-1);
mbacomdim_res.F.pvalall = 1 - fcdf(mbacomdim_res.F.valall, n-1, n-1);
mbacomdim_res.F.crit95 = finv(0.95, size(Y,1)-1, size(Y,1)-1);
mbacomdim_res.F.crit99 = finv(0.99, size(Y,1)-1, size(Y,1)-1);

% Note that the fields F.val/F.pval refer to the significance tests based
% on CC1 solely for the calculation of the F-ratios. 
% F.valall/F.pvalall refer to the use of all CCs linked to the residuals.

mbacomdim_res.Options = Options;

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    
    for i = 1 : length(Xs)
        tmp = mbacomdim_res.adecomp(i); 
        tmp.Y = mbacomdim_res.Y;
        tmp.Options = mbacomdim_res.Options;
        [mbacomdim_res.ptest(i)] = mfdptest(tmp);
    end
    
end


end