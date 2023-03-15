function [amopls_res] = amopls(X, Y, Options)
%
% The amopls function performs ANOVA – multiblock Orthogonal Partial
% Least Squares (AMOPLS). In this implementation, the ANOVA decomposition
% step can be peformed using different strategies. These strategies contain
% the classical ANOVA decomposition based on the differences of means, the
% ANOVA decomposition based on the general linear model methodology using
% sum or weighted-effect coding approaches. The last strategy is based on
% a rebalanced strategy. All these approaches give the same
% results when the experimental design is balanced but will differ when the
% design is unbalanced.
%
% AMOPLS makes use of the K-OPLS package for Matlab version 1.1.0
% https://sourceforge.net/projects/kopls/
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
%   'rebalanced';
%   Options.coding : if Options.decomp is 'glm', the coding method of the
%   design matrices can be chosen as 'sumcod' or 'wecod';
%   Options.interactions : the maximum number of interactions to consider
%   Options.permtest : permutation tests for effects significance if ~= 0;
%   Options.nperms : number of permutations to perform
%   Options.noc : number of orthogonal components to extract;
%
% Output arguments :
% ==================
% amopls_res : structure containing results of the AMOPLS
%   amopls_res.Xprepro : stores the preprocessed data before decomposition
%   amopls_res.Y: numerical array of the experimental design;
%   amopls_res.Y: numerical array of the experimental data;
%   amopls_res.Xm: the grand mean matrix;
%   amopls_res.Xf: cell array of the pure effect matrices;
%   amopls_res.Xe: numerical array of the pure error (residuals);
%   amopls_res.ssq: sum of squares of the effects;
%   amopls_res.ssqvarexp: sum of squares explained variation;
%   amopls_res.Xfaug: cell array of the augmented effect matrices;
%   amopls_res.effects: cell array of the effects (factor indices);
%
% If the ANOVA decomposition uses the GLM methodology, also contains:
%   amopls_res.Ef: the residuals without considering the effect f in the model;
%   amopls_res.cod: the coding of the experimental design Y;
%   amopls_res.B: the GLM parameters;
%
%   amopls_res.scores : global scores of AMOPLS
%   amopls_res.loadings : global loadings of AMOPLS
%   amopls_res.lambda : weight of each block in of AMOPLS
%   amopls_res.maxlambda : explained effect by each salience
%   amopls_res.varexp : explained variance of AMOPLS
%   amopls_res.kopls : All ComDim results, see comdimc.m for more information
%   amopls_res.Options : options used to perform AMOPLS
%
%   if Options.permtest == 1 :
%   amopls_res.ptest : stores permutation tests results and p-values if
%       performed
%
%   amopls_res.Options : options used to perform AMOPLS
%
% Usage :
% =======
% AMOPLS : classical ANOVA decomposition
% Options.decomp = 'classical';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [asca_res] = asca(X, Y, Options);
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
% Options.coding = 'wecod';
% Options.interactions = 2;
% Options.permtest = 1;
% Options.nperms = 1000;
% Options.noc = 2; % or 1:2 for noc optimization
% [amopls_res] = amopls(X, Y, Options);
%
% or
%
% RAMOPLS : Rebalanced AMOPLS ANOVA decomposition
% Options.decomp = 'rebalanced';
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
% radecomp.m (performs rebalanced ANOVA decomposition)
% sadecomp.m (performs synthetic ANOVA decomposition)
% mfdptest.m (performs permutation tests)
%
% KOPLS package functions :
%   - koplsModel.m
%   - koplsPredict.m
%
% Authors :
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
warning('off','all')

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

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : q
    [res] = classcol(Y(:,i));
    Y(:,i) = res.classvec;
end

%% Preprocessing of the input data

[amopls_res.Xprepro] = xpreproc(X, Options); % preprocessing as chosen

%% Performs ANOVA decomposition of each variable into effect matrices

% Decomposition is either classical, using GLM methodology or rebalancing procedure
if strcmp(Options.decomp,'rebalanced') == 1
    [adecomp_res] = radecomp(amopls_res.Xprepro.data, Y, Options); % see radecomp.m for information
else
    [adecomp_res] = adecomp(amopls_res.Xprepro.data, Y, Options); % see adecomp.m for information
end

% Copies fields of adecomp structure to the asca_res structure
fieldn = fieldnames(adecomp_res)';
for fn = fieldn(1:end-1) % without options
    amopls_res.(fn{1}) = adecomp_res.(fn{1});
end

%% Performs OPLS on the original design for each number of orthogonal components

% Calculates the association matrices
[Ws, Ycat] = association_matrices(adecomp_res);

% Calculates the kopls model based on sum of all association matrices for
% each number of orthogonal components defined in Options.noc
for i = 1 : length(Options.noc)
    if i > 1
        fieldn = fieldnames(amopls_res)';
        for fn = fieldn(1:13) % without options
            amopls_res(i).(fn{1}) = amopls_res(1).(fn{1});
        end
    end
    [amopls_res(i).kopls, amopls_res(i).lambda, ...
        amopls_res(i).scores, amopls_res(i).loadings, ...
        amopls_res(i).rsr,amopls_res(i).maxlambda] = kopls(adecomp_res, Ws, Ycat, Options.noc(i));
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
                    ptest.ssq(k,i,j) = perm_res.ssq(i);
                    ptest.r2y(k,i,j) = koplstmp.R2Yhat(end);
                    ptest.rsr(k,i,j) = rsr(i);

                else % permutation tests for the interaction terms

                    % Stores the permutation results
                    ptest.ssq(k,q+1:length(perm_res.effects),j) = perm_res.ssq(q+1:end-1);
                    ptest.r2y(k,q+1:length(perm_res.effects),j) = koplstmp.R2Yhat(end);
                    ptest.rsr(k,q+1:length(perm_res.effects),j) = rsr(q+1:end-1);

                end

            end
        end
    end

    % Calculates the p-values as the percentage of ssq based on
    % permutations equal or greater than the true experimental ssq
    % Performed similarly for rsr and r2y.
    for i = 1 : length(amopls_res(1).effects)
        for j = 1 : length(Options.noc) % loops over the number of orthogonal components
            ptest.pval.ssq(j,i) = sum(ptest.ssq(:,i,j) >= amopls_res(1).ssq(i)) / Options.nperms;
            ptest.pval.r2y(j,i) = sum(ptest.r2y(:,i,j) >= amopls_res(j).kopls.R2Yhat(end)) / Options.nperms;
            ptest.pval.rsr(j,i) = sum(ptest.rsr(:,i,j) >= amopls_res(j).rsr(i)) / Options.nperms;
        end
    end

    % Defines the optimal number of orthogonal components
    if Options.interactions == 1
        nOC = min(Options.noc(find(mean(ptest.pval.r2y(:,1:size(Y,2)),2) == min(mean(ptest.pval.r2y(:,1:size(Y,2)),2)))));
    else
        nOC = min(Options.noc(find(mean(ptest.pval.r2y(:,1:size(Y,2)+1),2) == min(mean(ptest.pval.r2y(:,1:size(Y,2)+1),2)))));
    end
    sfields = fieldnames(amopls_res);

    % Prepares the final structure of the AMOPLS results
    % Prepares the final structure of the AMOPLS results
    for i = 1 : length(sfields)
        amopls_res(length(Options.noc)+1).(sfields{i}) = amopls_res(nOC).(sfields{i});
    end

    amopls_res = amopls_res(length(Options.noc)+1);

    % Stores the nOCs to extract
    amopls_res.noc = nOC;

    % Stores the permutation test results
    amopls_res.ptest = ptest;

elseif length(Options.noc) == 1 % if only a scalar value was defined

    % Stores the nOCs to extract
    nOC = Options.noc;
    amopls_res.noc = nOC;

end

%% Calculates the percentage of explained variance of each predictive component

amopls_res.varexp = zeros(1, size(amopls_res.scores,2));
sstot_K = (sum(diag(amopls_res.kopls.K{1,1})));
for i = 1 : size(amopls_res.scores,2)
    rss = sum(diag( amopls_res.kopls.K{1,1} - (amopls_res.scores(:,i) * amopls_res.scores(:,i)') ));
    amopls_res.varexp(i) = (1 - rss/sstot_K) * 100;
end

%% In case Rebalanced AMOPLS was used, the algorithm enters this section

if strcmp(Options.decomp,'rebalanced') == 1

    % ANOVA decomposition of the original matrix according to the mean
    % effects of each factor and interaction
    [sadecomp_res] = sadecomp(amopls_res.Xprepro.data, Y, adecomp_res);

    % Copies fields of adecomp structure to the asca_res structure
    fieldn = fieldnames(sadecomp_res)';
    for fn = fieldn(1:end-1) % without options
        amopls_res.(fn{1}) = sadecomp_res.(fn{1});
    end

    % Augmented effect matrices + matrix of residuals
    Xs = [sadecomp_res.Xfaug,sadecomp_res.Xe]; % prediction
    Xsc = [adecomp_res.Xfaug,adecomp_res.Xe]; % calibration

    % Calculates the association matrices for the prediction of the scores
    for i = 1 : length(Xs)

        % Between the augmented effect matrices
        Wsp{i} = ( Xs{i} * Xs{i}' );
        Wsp{i} = Wsp{i} / norm(Wsp{i},'fro');

        % Between the augmented effect matrices and the Rebalanced
        % mean augmented effect matrices (needed by KOPLS package)
        Wtetr{i} = ( Xs{i} * Xsc{i}' );
        Wtetr{i} = Wtetr{i} / norm(Wtetr{i},'fro');

    end

    % Calculates the concatenated rebalanced scores and orthogonal scores
    [kopls_pred] = koplsPredict(sum(cat(3,Wtetr{:}),3),sum(cat(3,Wsp{:}),3),sum(cat(3,Ws{:}),3),amopls_res.kopls,nOC,0);
    amopls_res.scores = [kopls_pred.Tp{end},cat(2,kopls_pred.to{:})];
    amopls_res.kopls.adecomp = adecomp_res;
end

warning('on','all');

amopls_res.Options = Options; % stores options used

if strcmp(Options.decomp,'rebalanced') == 1
    amopls_res.Options.Xbal = adecomp_res.X;
end

%% Performs permutation tests if asked to in Options.permtest

% If the function performs permutation tests, it randomizes samples so that
% they are assigned to different levels in order to assert significance
if Options.permtest ~= 0 && length(Options.noc) == 1
    [amopls_res.ptest] = mfdptest(amopls_res);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to calculate AMOPLS association matrices
function [Ws, Ycat, neff] = association_matrices(adecomp_res)

% Number of observations
n = size(adecomp_res.X,1);

% Number of effects (including residuals)
neff = length(adecomp_res.effects);

% Allocates space for the association matrices
Ws = cell(1,neff);

% Allocates space for the kopls scores
Ycat = [];

% Loops over all effects + residuals
for i = 1 : neff+1

    % Enters the if statement if something else than residuals (last table)
    if i < neff+1 % effects

        [U,S,~] = svd(adecomp_res.Xf{i}); % SVD decomposition of each effect
        foo = diag(S); % diagonalizes singular values
        tol = max(size(adecomp_res.Xf{i})) * eps(norm(adecomp_res.Xf{i}));
        % r = find(foo < tol,1); % Finds non-zero LVs
        r = find(foo > tol,1,'last'); % Finds non-zero LVs
        % r = rank(adecomp_res.avgmat{i+1}) + 1;
        Ycat = [Ycat, U(:,1:r) .* repmat(foo(1:r)',n,1)];

        % Stores and normalizes to unit variance the association matrices
        Ws{i} = ( adecomp_res.Xfaug{i} * adecomp_res.Xfaug{i}' );

    else % residuals

        Ws{i} = ( adecomp_res.Xe * adecomp_res.Xe' );

    end

    Ws{i} = Ws{i} / norm(Ws{i},'fro');

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to perform KOPLS and side calculations
function [kopls, lambda, scores, loadings, rsr, maxsal] = kopls(adecomp_res, Ws, Ycat, nc)

Xs = [adecomp_res.Xfaug,adecomp_res.Xe];

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
    [~,maxsal(i)] = max(lambda(:,i));
end

% Concatenates scores and orthogonnal scores
scores = [kopls.T kopls.To];

% Allocates space for loadings
loadings = zeros(size(adecomp_res.X,2),size(scores,2));

% Stores loadings
for i = 1 : size(Ycat,2)
    loadings(:,i) = Xs{maxsal(i)}' * kopls.T(:,i) / (kopls.T(:,i)' * kopls.T(:,i));
end

% Stores orthogonal loadings
for i = 1 : nc
    loadings(:,size(Ycat,2)+i) = Xs{maxsal(size(Ycat,2)+i)}' * kopls.To(:,i) / (kopls.To(:,i)' * kopls.To(:,i));
end

% Calculation of the contribution ratio : noise to effect
for i = 1 : ntab
    rsr(i) = lambda(ntab,size(Ycat,2)+1)/lambda(i,size(Ycat,2)+1);
end

% Stores Y and raw lambda in kopls structure
kopls.Y = Ycat;
kopls.rawsaliences = lambda_raw;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inner function to perform permutations of specific effects.
% Used for significance testing and nOCs optimization.
function [adecomp_res] = amopls_perm(adecomp_res, effect, Options)

% Extracts the design matrices and the inital data
Y = adecomp_res.Y; % design matrix
X = adecomp_res.X; % experimental data

[n,p] = size(X); % size of X
[~,q] = size(Y); % size of Y

% Creates a copy of the design
Yvo = Y;

% Extracts the list of effects to include in the model
effects = adecomp_res.effects;
adecomp_res.effects = effects'; % stores the effect identifiers

% Switches cases depending on the multivariate ANOVA decomposition method
switch Options.decomp

    case 'classical'

        Y = Yvo; % extracts true design matrix

        if effect <= q % main effects

            % Reduces only to main effects
            effects = effects(1:q);
            adecomp_res.effects = effects;
            adecomp_res.Xf = adecomp_res.Xf(1:q);
            adecomp_res.Xfaug = adecomp_res.Xfaug(1:q);
            adecomp_res.ssq = adecomp_res.ssq(1:q);
            adecomp_res.ssqvarexp = adecomp_res.ssqvarexp(1:q);

            % Permutes only one main effect at a time if the effect input parameters
            % relates to one of the main factors
            Yo = Y(:,1:end ~= effects{effect}); % design matrix of other factors than effect
            Youni = unique(Yo,'rows'); % unique combinations of other factors

            % Creates randomized design matrix for factor effect without
            % affecting the other factors
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
                adecomp_res.Xf{1,effect}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
            end

            % Removes the grand mean from the effect matrix it is
            % calculated based on the initial X matrix
            adecomp_res.Xf{1,effect} = adecomp_res.Xf{1,effect} - adecomp_res.Xm;

        else % interactions

            % For interactions, the whole matrix is randomized
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
            Xred = adecomp_res.X - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{1,idx}),3);

            % Loops over combinations of levels in factors
            for j = 1 : length(uni)
                idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of observations of factors combs{i} at levels uni{j}
                levmean = mean(Xred(idx,:),1); % mean of observations of factors combs{i} at levels uni{j}
                adecomp_res.Xf{1,effect}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
            end

        end

        % Calculates the residuals
        adecomp_res.Xe = adecomp_res.X - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{:}),3);

    case 'glm'

        Y = Yvo; % extracts true design matrix

        X = adecomp_res.X; % experimental data

        if effect <= q % main effects

            Yo = Y(:,1:end ~= effect); % design matrix of other factors than i
            Youni = unique(Yo,'rows'); % unique combinations of other factors

            % Creates randomized design matrix for factor i without
            % affecting the other factors
            for j = 1 : size(Youni)
                idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                Y(idx,:) = Y(idx(perms),:); % performs the permutations here
            end

            % Rebalances the design
            [adecomp_res] = radecomp(X, Y, Options); % see radecomp.m for information

            % Reduces only to main effects
            effects = effects(1:q);
            adecomp_res.effects = effects;
            adecomp_res.Xf = adecomp_res.Xf(1:q);
            adecomp_res.Xfaug = adecomp_res.Xfaug(1:q);
            adecomp_res.ssq = adecomp_res.ssq(1:q);
            adecomp_res.ssqvarexp = adecomp_res.ssqvarexp(1:q);

        else

            % For interactions, the whole matrix is randomized
            perms = randperm(n); % random permutations
            Y = Yvo(perms,:); % randomizes levels assignment

            % Rebalances the design
            [adecomp_res] = radecomp(X, Y, Options); % see radecomp.m for information

        end

        Ys = adecomp_res.Y(:,adecomp_res.effects{effect});

        Ysm{1,1} = unique(Ys,'rows');

        for j = 1 : size(Ysm{1},1)

            [~,idx] = ismember(Ys,Ysm{1}(j,:), 'rows');
            Ysm{1,2}(j,:) = mean(adecomp_res.Xf{effect}(idx == 1, :),1);

        end

        for j = 1 : size(Ysm{1},1)
            idx = ismember(Y(:, unique(adecomp_res.effects{effect})), Ysm{1}(j,:),'rows');
            adecomp_res.Xf{effect}(idx,:) = repmat(Ysm{2}(j,:),[sum(idx==1),1]);
        end

        % Calculates the residuals
        adecomp_res.Xe = adecomp_res.X - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{:}),3);

    case 'rebalanced'

        Y = Yvo; % extracts true design matrix

        if effect <= q % main effects

            % Permutes only one main effect at a time if the effect input parameters
            % relates to one of the main factors
            Yo = Y(:,1:end ~= effects{effect}); % design matrix of other factors than effect
            Youni = unique(Yo,'rows'); % unique combinations of other factors

            % Creates randomized design matrix for factor effect without
            % affecting the other factors
            for j = 1 : size(Youni)
                idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                Y(idx,effects{effect}) = Y(idx(perms),effects{effect}); % randomizes the design for factor effect
            end

            % Rebalances the design
            [adecomp_res] = radecomp(X, Y, Options); % see radecomp.m for information

            % Reduces only to main effects
            effects = effects(1:q);
            adecomp_res.effects = effects;
            adecomp_res.Xf = adecomp_res.Xf(1:q);
            adecomp_res.Xfaug = adecomp_res.Xfaug(1:q);
            adecomp_res.ssq = adecomp_res.ssq(1:q);
            adecomp_res.ssqvarexp = adecomp_res.ssqvarexp(1:q);

        else % interactions

            % For interactions, the whole matrix is randomized
            perms = randperm(n); % random permutations
            Y = Yvo(perms,:); % randomizes levels assignment

            % Rebalances the design
            [adecomp_res] = radecomp(X, Y, Options); % see radecomp.m for information

        end

        Ys = adecomp_res.Y(:,adecomp_res.effects{effect});

        Ysm{1,1} = unique(Ys,'rows');

        for j = 1 : size(Ysm{1},1)

            [~,idx] = ismember(Ys,Ysm{1}(j,:), 'rows');
            Ysm{1,2}(j,:) = mean(adecomp_res.Xf{effect}(idx == 1, :),1);

        end

        for j = 1 : size(Ysm{1},1)
            idx = ismember(Y(:, unique(adecomp_res.effects{effect})), Ysm{1}(j,:),'rows');
            adecomp_res.Xf{effect}(idx,:) = repmat(Ysm{2}(j,:),[sum(idx==1),1]);
        end

        % Calculates the residuals
        adecomp_res.Xe = adecomp_res.X - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{:}),3);

end

% Calculates the sum of square and explained variation of the effects + residuals
neff = length(adecomp_res.Xf); % number of effects
for i = 1 : neff + 1

    if i < neff + 1

        % Sum-of-squares
        adecomp_res.ssq(1,i) = sum(adecomp_res.Xf{1,i}(:).^2);

        % Explained variation
        adecomp_res.ssqvarexp(1,i) = adecomp_res.ssq(1,i) / ( sum(adecomp_res.X(:).^2) - sum(adecomp_res.Xm(:).^2) ) * 100;

    else

        % Sum-of-squares
        adecomp_res.ssq(1,i) = sum(adecomp_res.Xe(:).^2);

        % Explained variation
        adecomp_res.ssqvarexp(1,i) = adecomp_res.ssq(1,i) / ( sum(adecomp_res.X(:).^2) - sum(adecomp_res.Xm(:).^2) ) * 100;

    end

end

% Calculates the augmented effect matrices Xfaug
adecomp_res.Xfaug = cell(1,length(adecomp_res.Xf));
for i = 1 : length(adecomp_res.Xf)
    adecomp_res.Xfaug{i} = adecomp_res.Xf{i} + adecomp_res.Xe;
end

% Adds the options field to the structure
adecomp_res.Options = Options;

end
