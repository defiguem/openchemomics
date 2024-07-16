function [iamopls_res] = iamopls(Xs, Yrep, Ys, Options)
%
% The amopls function performs ANOVA – multiblock Orthogonal Partial
% Least Squares (AMOPLS). In this implementation, the ANOVA decomposition
% step can be peformed using different strategies. These strategies contain
% the classical ANOVA decomposition based on the differences of means, the
% ANOVA decomposition based on the general linear model methodology using
% sum or weighted-effect coding approaches. The last strategy is based on
% a repeated undersampling strategy. All these approaches give the same
% results when the experimental design is balanced but will differ when the
% design is unbalanced.
%
% References :
% ============
% Boccard, J., & Rudaz, S. (2016). Exploring Omics data from designed
% experiments using analysis of variance multiblock Orthogonal Partial
% Least Squares. Analytica Chimica Acta, 920, 18‑28.
% https://doi.org/10.1016/j.aca.2016.03.042
%
% Bylesjö, M., Rantalainen, M., Nicholson, J. K., Holmes, E., & Trygg, J.
% (2008). K-OPLS package : Kernel-based orthogonal projections to latent
% structures for prediction and interpretation in feature space.
% BMC Bioinformatics, 9(1), 106. https://doi.org/10.1186/1471-2105-9-106
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
%   Options.noc : number of orthogonal components to extract;
%
% Output arguments :
% ==================
% amopls_res : structure containing results of the AMOPLS
%   amopls_res.Y : stores the design matrix
%   amopls_res.Xprepro : stores the preprocessed data before decomposition
%   amopls_res.anovamat : ANOVA matrices contained in cells (the last cell
%       contains the residuals)
%   amopls_res.avgmat : cell structure of factor/interation effects
%       (see amopls_res.matid for information on factor/interactions)
%   amopls_res.matid : identity of ANOVA matrices (grand mean, effects,
%       interactions and residual pure error)
%   amopls_res.ssq : sum of squares for each effect based on amopls_res.avgmat
%   amopls_res.ssqvarexp : variance explained by the sum of squares
%   amopls_res.Ys : stores the design matrices of factors and interactions
%   amopls_res.Xs : stores the effect matrix to which the pure error
%       was added
%   amopls_res.Ysm : mean effect vector for each factor and interaction
%   amopls_res.scores : AMOPLS scores for mean resampled ANOVA matrices
%   amopls_res.loadings : AMOPLS loadings for ANOVA matrices
%   amopls_res.pscores : projected AMOPLS scores for ANOVA matrices
%   amopls_res.Options : options used to perform AMOPLS
%
%   if Options.permtest == 1 :
%   amopls_res.ptest : stores permutation tests results and p-values if
%       performed
%
%   if Options.decomp == 'res_classical' || 'res_glm':
%   amopls_res.pavgmat : cell structure of factor/interation effects
%   amopls_res.panovamat : deflated ANOVA matrices contained in cells
%   amopls_res.pssq : sum of squares for each effect based on the synthetic
%   effect matrices
%
% Usage :
% =======
% AMOPLS : classical ANOVA decomposition
% Options.decomp = 'classical';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [amopls_res] = amopls(X, Y, Options);
%
% or
%
% AMOPLS+ : sum-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [amopls_res] = amopls(X, Y, Options);
%
% or
%
% WE-AMOPLS : weighted-effect-coding GLM ANOVA decomposition
% Options.decomp = 'glm';
% Options.coding = 'sumcod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [amopls_res] = amopls(X, Y, Options);
%
% or
%
% RAMOPLS : Resampled AMOPLS classical ANOVA decomposition
% Options.decomp = 'res_classical'; % or 'res_glm' but gives same results
% Options.nres = 100; % number of resamplings to perform
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [amopls_res] = amopls(X, Y, Options);
%
% Related function :
% ==================
% xpreproc.m (preprocessing of the input data matrix)
% adecomp.m (performs ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
% radecomp.m (performs multiple resamplings before ANOVA decomposition)
% sadecomp.m (ANOVA decomposition of the orignal matrix for Resampled AMOPLS)
% mfdptest.m (performs permutation tests)
%
% KOPLS package functions :
%   - koplsModel.m
%   - koplsPredict
%
% Authors :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Julien Boccard
% @ julien.boccard.unige.ch
%
% Modifications:
% ==============
%
% =========================================================================

%% Fail-safe section

warning('off','all')

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

% Checks if the number of orthogonal components (nOCs) to extract was
% defined. If not, it is set to one. If it is a scalar, it extracts the
% associated nOCs. If it contains more than one scalar then the nOCs need
% to be optimized optimized
if isfield(Options,'noc') == 0
    Options.noc = 1;
end

% Checks if permutation tests must be conducted
if length(Options.noc) > 1
    Options.permtest = 1;
elseif isfield(Options,'permtest') == 0
    Options.permtest = 0;
end

% Checks if the number of permutation tests to perform by effect was chosen
if isfield(Options,'nperms') == 0
    Options.nperms = 1000;
end

%% Preprocessing of the input data

[iamopls_res.Xprepro] = xpreproc(Xs, Options); % preprocessing as chosen

%% Performs multivariate ANOVA decomposition into effect matrices

% Note that for integrative approaches, rebalanced strategies should be
% favored unless all matrices have shared sample mode.
if strcmp(Options.decomp,'rebalanced') == 1
    [adecomp_res, Xs, Ys] = iradecomp({iamopls_res.Xprepro.data}, Yrep, Ys, Options); % see iradecomp.m for information
    Y = Ys{1};
else
    for i = 1 : length(Xs)
        [adecomp_res(i)] = adecomp(iamopls_res.Xprepro(i).data, Ys{1}, Options); % see adecomp.m for information
    end
    Y = Ys{1};
end

%% Performs OPLS on the original design for each number of orthogonal components

% Calculates the association matrices
[Ws, Ycat] = association_matrices(adecomp_res);

% Calculates the kopls model based on sum of all association matrices for
% each number of orthogonal components defined in Options.noc
for i = 1 : length(Options.noc)
    iamopls_res(i).Xprepro = iamopls_res(1).Xprepro;
    iamopls_res(i).X = Xs;
    iamopls_res(i).Y = Y;
    iamopls_res(i).adecomp = adecomp_res;
    [iamopls_res(i).kopls, iamopls_res(i).lambda, ...
        iamopls_res(i).scores, iamopls_res(i).loadings, ...
        iamopls_res(i).rsr,iamopls_res(i).maxlambda] = kopls(adecomp_res, Ws, Ycat, Options.noc(i));
end

%% Optimizes the number of orthogonal components if required

% If Options.noc contains more than one element, the optimization is
% performed to choose the optimal number of orthogonal components
if length(Options.noc) > 1

    for i = 1 : size(Y,2) + 1 % loops over the each main effect + once for all interactions
        for j = 1 : length(Options.noc) % loops over the number of orthogonal components
            for k = 1 : Options.nperms % loop over the number of permutations to perform

                clc; fprintf('Performing permutations test, please wait...\n Effect : %d / %d \n Orthogonal component : %d / %d \n Iteration : %d / %d \n', i, size(Y,2) + 1, j, length(Options.noc), k, Options.nperms);

                % Calculates the effect matrices based on permutations
                [perm_res] = amopls_perm(adecomp_res, i, Options);

                % Calculates the association matrices
                [Wstmp, Ycat] = association_matrices(perm_res);

                % Performs OPLS
                [koplstmp, ~, ~, ~,rsr] = kopls(perm_res, Wstmp, Ycat, Options.noc(j));

                if i <= q % permutation tests for the main effects

                    % Stores the permutation results
                    % ptest.ssq(k,i,j) = perm_res.ssq(i);
                    ptest.r2y(k,i,j) = koplstmp.R2Yhat(end);
                    % ptest.rsr(k,i,j) = rsr(i);

                else % permutation tests for the interaction terms

                    % Stores the permutation results
                    % ptest.ssq(k,q+1:length(perm_res.effects),j) = perm_res.ssq(q+1:end-1);
                    ptest.r2y(k,q+1,j) = koplstmp.R2Yhat(end); %ptest.r2y(k,q+1:length(perm_res.effects),j) = koplstmp.R2Yhat(end);
                    % ptest.rsr(k,q+1:length(perm_res.effects),j) = rsr(q+1:end-1);

                end

            end
        end
    end

    % Calculates the p-values as the percentage of ssq based on
    % permutations equal or greater than the true experimental ssq
    % Performed similarly for rsr and r2y.
    for i = 1 : (q+1)
        for j = 1 : length(Options.noc) % loops over the number of orthogonal components
            % ptest.pval.ssq(j,i) = sum(ptest.ssq(:,i,j) >= amopls_res(1).ssq(i)) / Options.nperms;
            ptest.pval.r2y(j,i) = sum(ptest.r2y(:,i,j) >= iamopls_res(j).kopls.R2Yhat(end)) / Options.nperms;
            % ptest.pval.rsr(j,i) = sum(ptest.rsr(:,i,j) >= amopls_res(j).rsr(i)) / Options.nperms;
        end
    end

    % Defines the optimal number of orthogonal components
    if Options.interactions == 1
        nOC = min(Options.noc(find(mean(ptest.pval.r2y(:,1:size(Y,2)),2) == min(mean(ptest.pval.r2y(:,1:size(Y,2)),2)))));
    else
        nOC = min(Options.noc(find(mean(ptest.pval.r2y(:,1:size(Y,2)+1),2) == min(mean(ptest.pval.r2y(:,1:size(Y,2)+1),2)))));
    end
    sfields = fieldnames(iamopls_res);

    % Prepares the final structure of the AMOPLS results
    for i = 1 : length(sfields)
        iamopls_res(length(Options.noc)+1).(sfields{i}) = iamopls_res(nOC).(sfields{i});
    end

    iamopls_res = iamopls_res(length(Options.noc)+1);

    % Stores the nOCs to extract
    iamopls_res.noc = nOC;

    % Stores the permutation test results
    iamopls_res.ptest = ptest;

elseif length(Options.noc) == 1 % if only a scalar value was defined

    % Stores the nOCs to extract
    nOC = Options.noc;
    iamopls_res.noc = nOC;

end

%% Calculates the percentage of explained variance of each predictive component

iamopls_res.varexp = zeros(1, size(iamopls_res.scores,2));
sstot_K = (sum(diag(iamopls_res.kopls.K{1,1})));
for i = 1 : size(iamopls_res.scores,2)
    rss = sum(diag( iamopls_res.kopls.K{1,1} - (iamopls_res.scores(:,i) * iamopls_res.scores(:,i)') ));    
    iamopls_res.varexp(i) = (1 - rss/sstot_K) * 100;
end

%% In case Rebalanced AMOPLS was used, the algorithm enters this section

if strcmp(Options.decomp,'rebalanced') == 1

    % ANOVA decomposition of the original matrix according to the mean
    % effects of each factor and interaction
    [sadecomp_res] = isadecomp(iamopls_res.X, Y, adecomp_res);
    iamopls_res.adecomp = sadecomp_res;
    [iamopls_res.adecomp.nans] = adecomp_res.nans;

    % Augmented effect matrices + matrix of residuals
    Xsall = {sadecomp_res.Xfaug}; 
    Xsallc = {adecomp_res.Xfaug}; 
    for i = 1 : length(Xsall)
        Xsall{i} = [Xsall{i},sadecomp_res(i).Xe]; 
        Xsallc{i} = [Xsallc{i},adecomp_res(i).Xe]; 
    end
    Xsall = cat(2,Xsall{:});
    Xsallc = cat(2,Xsallc{:});

    % Calculates the association matrices for the prediction of the scores
    for i = 1 : length(Xsall)

        % Between the augmented effect matrices
        Wsp{i} = ( Xsall{i} * Xsall{i}' );
        Wsp{i} = Wsp{i} / norm(Wsp{i},'fro');

        % Between the augmented effect matrices and the Rebalanced
        % mean augmented effect matrices (needed by KOPLS package)
        Wtetr{i} = ( Xsall{i} * Xsallc{i}' );
        Wtetr{i} = Wtetr{i} / norm(Wtetr{i},'fro');

    end

    % Calculates the concatenated rebalanced scores and orthogonal scores
    [kopls_pred] = koplsPredict(sum(cat(3,Wtetr{:}),3),sum(cat(3,Wsp{:}),3),sum(cat(3,Ws{:}),3),iamopls_res.kopls,nOC,0);
    iamopls_res.scores = [kopls_pred.Tp{end},cat(2,kopls_pred.to{:})];

end

warning('on','all');

iamopls_res.Options = Options; % stores options used

%% Reshapes the block weights for clarity

ntabs = length(adecomp_res(1).effects)+1; % effects + residuals
nblocks = length(Xs); % number of sources

% Lambda weights
iamopls_res.lambda = permute(reshape(iamopls_res.lambda,[ntabs,nblocks,size(iamopls_res.scores,2)]),[2,1,3]);

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0
    
    for i = 1 : length(Xs)
        tmp = iamopls_res.adecomp(i); 
        tmp.Y = iamopls_res.Y;
        tmp.Options = iamopls_res.Options;
        [iamopls_res.ptest2(i)] = mfdptest(tmp);
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to calculate AMOPLS association matrices
function [Ws, Ycat, neff] = association_matrices(adecomp_res)

% Number of observations
n = size(adecomp_res(1).Y,1);

% Number of effects (including residuals)
neff = length(adecomp_res(1).effects); 

% Allocates space for the kopls scores
Ycat = [];

Xf = {adecomp_res.Xf};
Xf = vertcat(Xf{:});

% Loops over all effects + residuals
for i = 1 : neff

    Xfcat = cat(2,Xf{:,i}); % concatenates all pure effect of all tables
    [U,S,~] = svd(Xfcat); % SVD decomposition of each effect
    foo = diag(S); % diagonalizes singular values
    tol = max(size(Xfcat)) * eps(norm(Xfcat));
    % r = find(foo < tol,1); % Finds non-zero LVs
    r = find(foo > tol,1,'last'); % Finds non-zero LVs
    % r = rank(adecomp_res.avgmat{i+1}) + 1;
    Ycat = [Ycat, U(:,1:r) .* repmat(foo(1:r)',n,1)];

end

% Allocates space for the association matrices
Ws = cell(1,neff);

Xf = {adecomp_res.Xfaug}; 
for i = 1 : length(Xf); Xf{i} = [Xf{i},adecomp_res(i).Xe]; end
Xf = cat(2,Xf{:});

% Calculates association matrices
for i = 1 : length(Xf)
        % Stores and normalizes to unit variance the association matrices
        Ws{i} = ( Xf{i} * Xf{i}' );
        Ws{i} = Ws{i} / norm(Ws{i},'fro');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to perform KOPLS and side calculations
function [kopls, lambda, scores, loadings, rsr, maxsal] = kopls(adecomp_res, Ws, Ycat, nc)

Xs = {adecomp_res.Xfaug}; 
for i = 1 : length(Xs); Xs{i} = [Xs{i},adecomp_res(i).Xe]; end
matid = repmat(1:length(adecomp_res(1).effects)+1,[1,length(Xs)]);
Xs = cat(2,Xs{:});

ntab = length(Xs);

% Calibrates K-OPLS model by summing all associate matrices
kopls = koplsModel(sum(cat(3,Ws{:}),3), Ycat, size(Ycat,2), nc, 'mc', 'mc');

% Calculates the saliences for each block
lambda = zeros(ntab,size(Ycat,2));
for j = 1 : ntab

    % For the scores
    for k = 1 : size(Ycat,2)
        lambda(j,k) = kopls.T(:,k)' * Ws{j} * kopls.T(:,k);
    end

    % For the orthogonal scores
    for l = 1 : nc
        lambda(j,k+l) = kopls.To(:,l)' * Ws{j} * kopls.To(:,l);
    end

end

% Creates a copy of the non-normalized salience values
lambda_raw = lambda;

% Normalizes the saliences and finds the block with maximum salience
for i = 1 : size(lambda,2)
    lambda(:,i) = lambda(:,i)/sum(lambda(:,i));
    [~,id] = max(lambda(:,i));
    maxsal(i) = matid(id);
end

% Concatenates scores and orthogonnal scores
scores = [kopls.T kopls.To];

% Allocates space for loadings
loadings = cell(size(adecomp_res,2),size(Ycat,2));

Xs = {adecomp_res.Xfaug}; 
for i = 1 : length(Xs); Xs{i} = [Xs{i},adecomp_res(i).Xe]; end
Xs = vertcat(Xs{:});

% Stores loadings
for i = 1 : size(Ycat,2)
    for j = 1 : size(adecomp_res,2)
        loadings{j,i} = Xs{j,maxsal(i)}' * kopls.T(:,i) / (kopls.T(:,i)' * kopls.T(:,i));
    end
end

% Stores orthogonal loadings
for i = 1 : nc
    for j = 1 : size(adecomp_res,2)
        loadings{j,size(Ycat,2)+i} = Xs{j,maxsal(size(Ycat,2)+i)}' * kopls.To(:,i) / (kopls.To(:,i)' * kopls.To(:,i));
    end
end

% Calculation of the contribution ratio : noise to effect
lam = reshape(lambda(:,size(Ycat,2)+1),[max(matid),size(adecomp_res,2)])';
rsr = repmat(lam(:,end),[1,size(lam,2)]) ./ lam;

% Stores Y and raw lambda in kopls structure
kopls.Y = Ycat;
kopls.rawsaliences = lambda_raw;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to perform permutations of specific effects.
% Used for significance testing and nOCs optimization.
function [adecomp_res] = amopls_perm(adecomp_res, effect, Options)

for t = 1 : length(adecomp_res)

    clear Ys Ysm
    
    % Extracts the design matrices and the inital data
    Y = adecomp_res(t).Y; % design matrix
    X = adecomp_res(t).X; % experimental data

    [n,p] = size(X); % size of X
    [~,q] = size(Y); % size of Y

    % Creates a copy of the design
    Yvo = Y;

    % Extracts the list of effects to include in the model
    effects = adecomp_res(t).effects;
    adecomp_res(t).effects = effects'; % stores the effect identifiers

    % Switches cases depending on the multivariate ANOVA decomposition method
    switch Options.decomp

        case 'classical'

            Y = Yvo; % extracts true design matrix

            if effect <= q % main effects

                % Reduces only to main effects
                effects = effects(1:q);
                adecomp_res(t).effects = effects;
                adecomp_res(t).Xf = adecomp_res(t).Xf(1:q);
                adecomp_res(t).Xfaug = adecomp_res(t).Xfaug(1:q);
                adecomp_res(t).ssq = adecomp_res(t).ssq(1:q);
                adecomp_res(t).ssqvarexp = adecomp_res(t).ssqvarexp(1:q);

                % Permutes only one main effect at a time if the effect input parameters
                % relates to one of the main factors
                Yo = Y(:,1:end ~= effects{effect}); % design matrix of other factors than effect
                Youni = unique(Yo,'rows'); % unique combinations of other factors

                % Creates randomized design matrix for factor effect without
                % affecting the other factors
                rng(12345);
                for j = 1 : size(Youni)
                    idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                    perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                    Y(idx,effects{effect}) = Y(idx(perms),effects{effect}); % randomizes the design for factor effect
                end

                % Stores unique levels for factor Y(:,i)
                levels = unique(Y(:,effects{effect}));

                % Loops over the levels of a given factor
                for j = 1 : length(levels)
                    idx = find(Y(:,effects{effect}) == levels(j)); % indices of observations at factor i and level j
                    levmean = mean(X(idx,:),1); % mean of observations of factor i at level j
                    adecomp_res(t).Xf{1,effect}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                end

                % Removes the grand mean from the effect matrix it is
                % calculated based on the initial X matrix
                adecomp_res(t).Xf{1,effect} = adecomp_res(t).Xf{1,effect} - adecomp_res(t).Xm;

            else % interactions

                % For interactions, the whole matrix is randomized
                rng(12345);
                perms = randperm(n); % random permutations
                Y = Yvo(perms,:); % randomizes levels assignment

                % Combination of factors (interactions)
                Ycomb = cellstr(num2str(Y(:,effects{effect})));

                % Unique combinations of levels in factors
                uni = unique(Ycomb);

                % Identifies all factors and interactions involved in the
                % evaluation of the interaction combs{i} that are  lower
                % than length of combs{i}
                combs = {};
                for j = 1 : length(effects{effect})-1
                    combs = [combs; num2cell(nchoosek(effects{effect}, j),2)];
                end

                % Identifies the position of the above effects in order to
                % remove them below
                idx = [];
                for j = 1 : size(combs,1)
                    idx = [idx; find(cell2mat(cellfun(@(x) isequal(x,combs{j}),effects,'UniformOutput',false)) == 1)];
                end

                % Removes all the effect matrices identified in combs
                % according to their position defined in idx.
                % Note that the interaction term is calculated based on the
                % initial matrix X, then we remove the above identified
                % effect matrices + the grand mean.
                Xred = adecomp_res(t).X - adecomp_res(t).Xm - sum(cat(3,adecomp_res(t).Xf{1,idx}),3);

                % Loops over combinations of levels in factors
                for j = 1 : length(uni)
                    idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of observations of factors combs{i} at levels uni{j}
                    levmean = mean(Xred(idx,:),1); % mean of observations of factors combs{i} at levels uni{j}
                    adecomp_res(t).Xf{1,effect}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                end

            end

            % Calculates the residuals
            adecomp_res(t).Xe = adecomp_res(t).X - adecomp_res(t).Xm - sum(cat(3,adecomp_res(t).Xf{:}),3);

        case 'glm'

            Y = Yvo; % extracts true design matrix

            X = adecomp_res(t).X; % experimental data

            if effect <= q % main effects

                Yo = Y(:,1:end ~= effect); % design matrix of other factors than i
                Youni = unique(Yo,'rows'); % unique combinations of other factors

                % Creates randomized design matrix for factor i without
                % affecting the other factors
                rng(12345);
                for j = 1 : size(Youni)
                    idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                    perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                    Y(idx,:) = Y(idx(perms),:); % performs the permutations here
                end

                % Rebalances the design
                [adecomp_res(t)] = radecomp(X, Y, Options); % see radecomp.m for information

                % Reduces only to main effects
                effects = effects(1:q);
                adecomp_res(t).effects = effects;
                adecomp_res(t).Xf = adecomp_res(t).Xf(1:q);
                adecomp_res(t).Xfaug = adecomp_res(t).Xfaug(1:q);
                adecomp_res(t).ssq = adecomp_res(t).ssq(1:q);
                adecomp_res(t).ssqvarexp = adecomp_res(t).ssqvarexp(1:q);

            else

                % For interactions, the whole matrix is randomized
                rng(12345);
                perms = randperm(n); % random permutations
                Y = Yvo(perms,:); % randomizes levels assignment

                % Rebalances the design
                [adecomp_res(t)] = radecomp(X, Y, Options); % see radecomp.m for information

            end

            Ys = adecomp_res(t).Y(:,adecomp_res(t).effects{effect});

            Ysm{1,1} = unique(Ys,'rows');

            for j = 1 : size(Ysm{1},1)

                [~,idx] = ismember(Ys,Ysm{1}(j,:), 'rows');
                Ysm{1,2}(j,:) = mean(adecomp_res(t).Xf{effect}(idx == 1, :),1);

            end

            for j = 1 : size(Ysm{1},1)
                idx = ismember(Y(:, unique(adecomp_res(t).effects{effect})), Ysm{1}(j,:),'rows');
                adecomp_res(t).Xf{effect}(idx,:) = repmat(Ysm{2}(j,:),[sum(idx==1),1]);
            end

            % Calculates the residuals
            adecomp_res(t).Xe = adecomp_res(t).X - adecomp_res(t).Xm - sum(cat(3,adecomp_res(t).Xf{:}),3);

        case 'rebalanced'

            Y = Yvo; % extracts true design matrix

            if effect <= q % main effects

                % Permutes only one main effect at a time if the effect input parameters
                % relates to one of the main factors
                Yo = Y(:,1:end ~= effects{effect}); % design matrix of other factors than effect
                Youni = unique(Yo,'rows'); % unique combinations of other factors

                % Creates randomized design matrix for factor effect without
                % affecting the other factors
                rng(12345);
                for j = 1 : size(Youni)
                    idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                    perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                    Y(idx,effects{effect}) = Y(idx(perms),effects{effect}); % randomizes the design for factor effect
                end

                % Rebalances the design
                tmp = radecomp(X, Y, Options); % see radecomp.m for information
                for fn = fieldnames(tmp)'
                    adecomp_res(t).(fn{1}) = tmp.(fn{1});
                end

                % Reduces only to main effects
                effects = effects(1:q);
                adecomp_res(t).effects = effects;
                adecomp_res(t).Xf = adecomp_res(t).Xf(1:q);
                adecomp_res(t).Xfaug = adecomp_res(t).Xfaug(1:q);
                adecomp_res(t).ssq = adecomp_res(t).ssq(1:q);
                adecomp_res(t).ssqvarexp = adecomp_res(t).ssqvarexp(1:q);

            else % interactions

                % For interactions, the whole matrix is randomized
                rng(12345);
                perms = randperm(n); % random permutations
                Y = Yvo(perms,:); % randomizes levels assignment

                % Rebalances the design
                tmp = radecomp(X, Y, Options); % see radecomp.m for information
                for fn = fieldnames(tmp)'
                    adecomp_res(t).(fn{1}) = tmp.(fn{1});
                end

            end

            Ys = adecomp_res(t).Y(:,adecomp_res(t).effects{effect});

            Ysm{1,1} = unique(Ys,'rows');

            for j = 1 : size(Ysm{1},1)

                [~,idx] = ismember(Ys,Ysm{1}(j,:), 'rows');
                Ysm{1,2}(j,:) = mean(adecomp_res(t).Xf{effect}(idx == 1, :),1);

            end

            for j = 1 : size(Ysm{1},1)
                idx = ismember(Y(:, unique(adecomp_res(t).effects{effect})), Ysm{1}(j,:),'rows');
                adecomp_res(t).Xf{effect}(idx,:) = repmat(Ysm{2}(j,:),[sum(idx==1),1]);
            end

            % Calculates the residuals
            adecomp_res(t).Xe = adecomp_res(t).X - adecomp_res(t).Xm - sum(cat(3,adecomp_res(t).Xf{:}),3);

    end

    % Calculates the sum of square and explained variation of the effects + residuals
    neff = length(adecomp_res(t).Xf); % number of effects
    for i = 1 : neff + 1

        if i < neff + 1

            % Sum-of-squares
            adecomp_res(t).ssq(1,i) = sum(adecomp_res(t).Xf{1,i}(:).^2);

            % Explained variation
            adecomp_res(t).ssqvarexp(1,i) = adecomp_res(t).ssq(1,i) / ( sum(adecomp_res(t).X(:).^2) - sum(adecomp_res(t).Xm(:).^2) ) * 100;

        else

            % Sum-of-squares
            adecomp_res(t).ssq(1,i) = sum(adecomp_res(t).Xe(:).^2);

            % Explained variation
            adecomp_res(t).ssqvarexp(1,i) = adecomp_res(t).ssq(1,i) / ( sum(adecomp_res(t).X(:).^2) - sum(adecomp_res(t).Xm(:).^2) ) * 100;

        end

    end

    % Calculates the augmented effect matrices Xfaug
    adecomp_res(t).Xfaug = cell(1,length(adecomp_res(t).Xf));
    for i = 1 : length(adecomp_res(t).Xf)
        adecomp_res(t).Xfaug{i} = adecomp_res(t).Xf{i} + adecomp_res(t).Xe;
    end

    % Adds the options field to the structure
    adecomp_res(t).Options = Options;

end

end

