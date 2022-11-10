function [adecomp_res] = adecomp(X, Y, Options)
%
% The function adecomp performs the multivariate ANalysis Of VAriance
% (ANOVA) decomposition of the X matrix. X contains experiments in the rows
% and variables in the columns. It is usually the result of designed and
% organized experiments aiming to tackle specific research questions.
%
% To do so, a set of independent variables (factors) which are expected to
% contribute significantly to the experimental variability is chosen and
% levels of these factors are defined in order to perform experiments under
% controlled conditions. The results of the experiments, usually
% multivariate, are investigated to assert whether such factors have an
% effect on the system under study.
%
% The decomposition can be done in two ways :
% (1) classical ANOVA decomposition into effect matrices based on
% differences of means.
% (2) ANOVA decomposition based on the General Linear Models (GLM).
%
% References :
% ============
% Thiel, M., Féraud, B., & Govaerts, B. (2017). ASCA+ and APCA+:
% Extensions of ASCA and APCA in the analysis of unbalanced multifactorial
% designs: Analyzing unbalanced multifactorial designs with ASCA+ and APCA+.
% Journal of Chemometrics, 31(6), e2895. https://doi.org/10.1002/cem.2895
%
% Ali, N., Jansen, J., van den Doel, A., Tinnevelt, G. H., & Bocklitz, T.
% (2020). WE-ASCA: The Weighted-Effect ASCA for Analyzing Unbalanced
% Multifactorial Designs?A Raman Spectra-Based Example. Molecules, 26(1).
% https://doi.org/10.3390/molecules26010066
%
% Nieuwenhuis, R., Grotenhuis, M. te, & Pelzer, B. (2017). Weighted Effect
% Coding for Observational Data with wec. The R Journal, 9(1), 477?485.
%
% Input arguments :
% =================
% X : data matrix with observations in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.decomp : ANOVA decomposition method 'classical', 'glm',
%       'rebalanced'.
%   Options.coding : if Options.decomp is 'glm', the coding
%       method of the design matrices can be chosen as 'sumcod' or 'wecod'.
%   Options.interactions : the maximum number of interactions to consider
%
% Output arguments :
% ==================
% adecomp_res : structure containing results of the ANOVA decomposition
%   adecomp_res.Y: numerical array of the experimental design;
%   adecomp_res.Y: numerical array of the experimental data;
%   adecomp_res.Xm: the grand mean matrix;
%   adecomp_res.Xf: cell array of the pure effect matrices;
%   adecomp_res.Xe: numerical array of the pure error (residuals);
%   adecomp_res.ssq: sum of squares of the effects;
%   adecomp_res.ssqvarexp: sum of squares explained variation;
%   adecomp_res.Xfaug: cell array of the augmented effect matrices;
%   adecomp_res.effects: cell array of the effects (factor indices);
%   adecomp_res.Options: the options used for calculation;
%
% If the ANOVA decomposition uses the GLM methodology, also contains:
%   adecomp_res.Ef: the residuals without considering the effect f in the model;
%   adecomp_res.cod: the coding of the experimental design Y;
%   adecomp_res.B: the GLM parameters;
%
% Usage :
% =======
% Options.decomp = 'classical';
% Options.interactions = 2;
% [adecomp_res] = adecomp(X, Y, Options)
%
% or
%
% Options.decomp = 'glm';
% Options.coding = 'sumcod'; % or 'wecod'
% Options.interactions = 2;
% [adecomp_res] = adecomp(X, Y, Options)
%
% Related functions :
% ===================
% radecomp.m (performs rebalanced ANOVA decomposition)
% sadecomp.m (performs synthetic ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
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

% Checks if the maximum number of interactions to consider was defined
if isfield(Options,'interactions') == 0
    Options.interactions = q;
elseif isfield(Options,'interactions') == 1 && Options.interactions > 3
    Options.interactions = q; % maximum number of interactions cannot be > q
end

% Checks if the ANOVA decomposition method was specified
if isfield(Options,'decomp') == 0
    Options.decomp = 'glm';
end

% Checks if the ANOVA decomposition was set to 'rebalanced'
if isfield(Options,'decomp') == 1 && strcmp(Options.decomp,'rebalanced') == 1 && isfield(Options,'coding') == 0
    Options.coding = 'sumcod';
end

% Checks if the ANOVA decomposition was set to GLM, then checks if the
% coding for the design matrices was chosen
if isfield(Options,'decomp') == 1 && strcmp(Options.decomp,'glm') == 1 && isfield(Options,'coding') == 0
    Options.coding = 'sumcod';
end

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : q
    Y(:,i) = yformatconv(Y(:,i),'intvec');
end

%% Multivariate ANOVA decomposition

% Stores the design matrices and the inital data
adecomp_res.Y = Y; % design matrix
adecomp_res.X = X; % experimental data

% Allocates fields in the output structure
adecomp_res.Xm = {}; % grand mean matrix
adecomp_res.Xf = {}; % pure effect matrices
adecomp_res.Xe = {}; % error matrix of residuals
adecomp_res.Ef = {}; % residuals of a model without effect f
adecomp_res.ssq = []; % effect sum of squares
adecomp_res.ssqvarexp = []; % explained variation of sum of squares
adecomp_res.Xfaug = {}; % augmented effect matrices

% Defines the list of effects to include in the model
effects = {};
for i = 1 : Options.interactions
    effects = [effects; num2cell(nchoosek(1:q, i),2)];
end
adecomp_res.effects = effects'; % stores the effect identifiers

% Switches cases depending on the multivariate ANOVA decomposition method
switch Options.decomp

    case 'classical'

        adecomp_res.Xm = repmat(mean(X,1), n,1); % grand mean matrix

        % Loops over the effects
        for i = 1 : length(effects)

            if length(effects{i}) == 1 % main effects

                % Stores unique levels for factor Y(:,i)
                levels = unique(Y(:,effects{i}));

                % Calculates pure effects for each level j of factor i
                adecomp_res.Xf{1,i} = zeros(n,p);

                % Loops over the levels of a given factor
                for j = 1 : length(levels)
                    idx = find(Y(:,effects{i}) == levels(j)); % indices of observations at factor i and level j
                    levmean = mean(X(idx,:),1); % mean of observations of factor i at level j
                    adecomp_res.Xf{1,i}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                end

                % Removes the grand mean from the effect matrix it is
                % calculated based on the initial X matrix
                adecomp_res.Xf{1,i} = adecomp_res.Xf{1,i} - adecomp_res.Xm;

            else % interactions

                % Combination of factors (interactions)
                Ycomb = cellstr(num2str(Y(:,effects{i})));

                % Unique combinations of levels in factors
                uni = unique(Ycomb);

                % Calculates pure effects for each combination of levels
                adecomp_res.Xf{1,i} = zeros(n,p);

                % Loops over combinations of levels in factors
                for j = 1 : length(uni)
                    idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of observations of factors combs{i} at levels uni{j}
                    levmean = mean(X(idx,:),1); % mean of observations of factors combs{i} at levels uni{j}
                    adecomp_res.Xf{1,i}(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                end

                % Identifies all factors and interactions involved in the
                % evaluation of the interaction combs{i} that are  lower
                % than length of combs{i}
                combs = {};
                for j = 1 : length(effects{i})-1
                    combs = [combs; num2cell(nchoosek(effects{i}, j),2)];
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
                adecomp_res.Xf{1,i} = adecomp_res.Xf{1,i} - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{1,idx}),3);

            end
        end

        % Calculates the residuals
        adecomp_res.Xe = adecomp_res.X - adecomp_res.Xm - sum(cat(3,adecomp_res.Xf{:}),3);

        % Calculates the sum of square and explained variation of the effects + residuals
        neff = length(effects); % number of effects
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

        % Removes the Ef field because not used in classical adecomp
        adecomp_res = rmfield(adecomp_res,'Ef');

    case {'glm','rebalanced'}

        % Checks which coding is to be used to define the design matrices
        if strcmp(Options.coding,'sumcod') == 1

            % Sum coding of the input design matrix according to
            % Thiel et al. (2017)
            [cod, codempty] = sumcoding(Y, Options);
            adecomp_res.cod = cod;

        elseif strcmp(Options.coding,'wecod') == 1

            % Weighted-effect coding of the input design matrix according
            % to Nieuwenhuis et al. (2017)
            [cod, codempty] = wecoding(Y, Options);
            adecomp_res.cod = cod;

        end

        % Calculates the model parameters in B
        codcat = cat(2,cod{:});
        Bmat = (pinv(codcat' * codcat) * codcat' * X)';

        % Creates empty B and separates Bmat in array of cells (it is just useful for later calculations...)
        bempty = {};
        B = {};
        for i = 1 : length(cod)
            B{i} = Bmat(:,size(cat(2,cod{1:i-1}),2) + 1 : size(cat(2,cod{1:i}),2));
            bempty{i} = zeros(p,size(cod{i},2));
        end

        % Stores the parameters B
        adecomp_res.B = B;

        % Calculates full model Residuals
        adecomp_res.Xe = X - (cat(2,cod{:}) * Bmat');

        % Calculation of the main factor effects :
        for i = 1 : length(cod) % all main effects and interactions

            % Calculates the effect matrices
            Xf = codempty; Xf{i} = cod{i};
            Bf = bempty; Bf{i} = B{i};

            if i == 1 % grand mean matrix

                adecomp_res.Xm = cat(2,Xf{:}) * cat(2,Bf{:})';
            else
                 % effects

                % Stores pure effect matrix
                adecomp_res.Xf{1,i-1} = cat(2,Xf{:}) * cat(2,Bf{:})';


                % Calculates the residuals of the model with effect i
                Xf = cod; Xf{i} = zeros(size(cod{i}));
                Bf = B; Bf{i} = zeros(size(B{i}));
                adecomp_res.Ef{1,i-1} = X - (cat(2,Xf{:}) * cat(2,Bf{:})');

            end

        end

        % Sum of squares and explained variance 
        neff = length(cod)-1;
        den = sum(cat(3,adecomp_res.Ef{1:neff}).^2,'all') - (neff-1) * sum(adecomp_res.Xe(:).^2);
        for i = 1 : neff + 1

            if i < neff + 1 % effects

                % Sum of squares
                adecomp_res.ssq(1,i) = sum(adecomp_res.Ef{i}(:).^2) - sum(adecomp_res.Xe(:).^2);

                % Explained variance (row 1 does not sum to 100%, row 2 sums to 100%)
                adecomp_res.ssqvarexp(1,i) = adecomp_res.ssq(1,i) / ( sum((X-adecomp_res.Xm).^2,'all') ) * 100;
                adecomp_res.ssqvarexp(2,i) = adecomp_res.ssq(1,i) / ( den ) * 100;

            else % residuals

                % Sum of squares
                adecomp_res.ssq(1,i) = sum(adecomp_res.Xe(:).^2);

                % Explained variance (row 1 does not sum to 100%, row 2 sums to 100%)
                adecomp_res.ssqvarexp(1,i) = sum(adecomp_res.Xe(:).^2) / ( sum((X-adecomp_res.Xm).^2,'all') ) * 100;
                adecomp_res.ssqvarexp(2,i) = sum(adecomp_res.Xe(:).^2) / ( den ) * 100;

            end

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

