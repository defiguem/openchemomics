function [ptest_res] = mfdptest(mfda_res)
%
% The mfdptest function performs permutation tests in order to
% determine main effects/interactions significance. To do so, repeated
% randomizations of the rows in the design matrix are performed before
% calculating effect matrices.
%
% Whereas for interaction the rows of the design matrix are permuted
% without restrictions, for the main effects there is one. Permutations
% for one factor are only allowed within the levels of other factors. For
% more details on permutation tests, see Anderson et Braak, 2003.
%
% References :
% ============
% Vis, D. J., Westerhuis, J. A., Smilde, A. K., & van der Greef, J. (2007).
% Statistical validation of megavariate effects in ASCA.
% BMC Bioinformatics, 8(1). https://doi.org/10.1186/1471-2105-8-322
%
% Anderson, M., & Braak, C. T. (2003). Permutation tests for
% multi-factorial analysis of variance. Journal of Statistical Computation
% and Simulation, 73(2), 85-113. https://doi.org/10.1080/00949650215733
%
% Thiel, M., Féraud, B., & Govaerts, B. (2017). ASCA+ and APCA+ :
% Extensions of ASCA and APCA in the analysis of unbalanced multifactorial
% designs: Analyzing unbalanced multifactorial designs with ASCA+ and APCA+.
% Journal of Chemometrics, 31(6), e2895. https://doi.org/10.1002/cem.2895
%
% te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R.,
% Schmidt-Catran, A., & Konig, R. (2017). When size matters : Advantages of
% weighted effect coding in observational studies. International Journal of
% Public Health, 62(1), 163‑167. https://doi.org/10.1007/s00038-016-0901-1
%
% de Figueiredo, M., Giannoukos, S., Rudaz, S., Zenobi, R., & Boccard, J. 
% (s. d.). Efficiently handling high-dimensional data from multifactorial 
% designs with unequal group sizes using Rebalanced ASCA (RASCA). 
% Journal of Chemometrics, n/a(n/a), e3401. https://doi.org/10.1002/cem.3401
%
% Input arguments :
% =================
% mfda_res : output argument of the adecomp.m function or of any other
% supported function performing analysis of multifactorial designs.
%
% Output arguments :
% ==================
% ptest_res : structure containing results of the permutation test
%   ptest_res.ssq : sum of squares from permutation tests (null-hypothesis)
%   ptest_res.tssq : experimental sum of squares
%   ptest_res.pval : p-values for each factor and interaction
%
% Usage :
% =======
% [ptest_res] = mfdptest(mfda_res);
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

% Retrieves options field from mfda_res structure
Options = mfda_res.Options;

% Checks if the number of permutation tests to perform by effect was chosen
if isfield(Options,'nperms') == 0
    Options.nperms = 1000;
end

if strcmp(Options.decomp,'glm') == 1 && isfield(Options,'coding') == 0
    Options.coding = 'sumcod';
end

%% Permutation tests

X = mfda_res.X; % experimental data
[n,p] = size(X); % size of X
[~,q] = size(mfda_res.Y); % size of Y

% Allocates space for ssq of permuted levels
ptest_res.ssq = zeros(Options.nperms, length(mfda_res.effects));

% Defines rng seed for repeatable results
rng(1000);

% Extracts the list of effects
effects = mfda_res.effects;

switch Options.decomp

    case 'classical' % classical ANOVA decomposition

        for k = 1 : Options.nperms

            clc; fprintf('Performing permutations test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);

            Xm = repmat(mean(X,1), n,1); % grand mean matrix

            % Loops over the effects
            for i = 1 : length(effects)

                if length(effects{i}) == 1 % main effects

                    Y = mfda_res.Y; % extracts true design matrix
                    Yo = Y(:,1:end ~= i); % design matrix of other factors than i
                    Youni = unique(Yo,'rows'); % unique combinations of other factors

                    % Creates randomized design matrix for factor i without
                    % affecting the other factors
                    for j = 1 : size(Youni)
                        idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                        perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                        Y(idx,i) = Y(idx(perms),i); % randomizes the design for factor i
                    end

                    % Stores unique levels for factor Y(:,i)
                    levels = unique(Y(:,effects{i}));

                    % Loops over the levels of a given factor
                    Xf = zeros(size(X));
                    for j = 1 : length(levels)
                        idx = find(Y(:,effects{i}) == levels(j)); % indices of observations at factor i and level j
                        levmean = mean(X(idx,:),1); % mean of observations of factor i at level j
                        Xf(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                    end

                    % Removes the grand mean from the effect matrix 
                    Xf = Xf - Xm;

                    % Calculates the sum of squares
                    ptest_res.ssq(k,i) = sum(Xf(:).^2);

                else % interactions

                    % For interactions, the whole matrix is randomized
                    perms = randperm(n); % random permutations
                    Y = mfda_res.Y(perms,:); % randomizes levels assignment

                    % Combination of factors (interactions)
                    Ycomb = cellstr(num2str(Y(:,effects{i})));

                    % Unique combinations of levels in factors
                    uni = unique(Ycomb);

                    % Loops over combinations of levels in factors
                    % For the interactions, we remove first the grand meand
                    % and all mains effects.
                    Xred = X - Xm - sum(cat(3,mfda_res.Xf{1:q}),3);

                    % Loops over the levels of a given interaction
                    Xf = zeros(size(X));
                    for j = 1 : length(uni)
                        idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of observations of factors combs{i} at levels uni{j}
                        levmean = mean(Xred(idx,:),1); % mean of observations of factors combs{i} at levels uni{j}
                        Xf(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                    end

                    % Calculates the sum of squares
                    ptest_res.ssq(k,i) = sum(Xf(:).^2);

                end
            end
        end

    case 'glm'

        for k = 1 : Options.nperms

            clc; fprintf('Performing permutations test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);

            % Loops over the effects
            for i = 1 : length(effects)

                if length(effects{i}) == 1 % main effects

                    Y = mfda_res.Y; % extracts true design matrix
                    Yo = Y(:,1:end ~= i); % design matrix of other factors than i
                    Youni = unique(Yo,'rows'); % unique combinations of other factors

                    % Creates randomized design matrix for factor i without
                    % affecting the other factors
                    for j = 1 : size(Youni)
                        idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                        perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                        Y(idx,i) = Y(idx(perms),i); % randomizes the design for factor i
                    end

                    % Checks which coding is to be used to define the design matrices
                    if strcmp(Options.coding,'sumcod') == 1

                        % Sum coding of the input design matrix according to
                        % Thiel et al. (2017)
                        [cod, codempty] = sumcoding(Y, Options);

                    elseif strcmp(Options.coding,'wecod') == 1

                        % Weighted-effect coding of the input design matrix according
                        % to Nieuwenhuis et al. (2017)
                        [cod, codempty] = wecoding(Y, Options);

                    end

                    % Calculates the model parameters in B
                    codcat = cat(2,cod{:});
                    Bmat = (pinv(codcat' * codcat) * codcat' * X)';

                    % Creates empty B and separates Bmat in array of cells (it is just useful for later calculations...)
                    bempty = {};
                    B = {};
                    for j = 1 : length(cod)
                        B{j} = Bmat(:,size(cat(2,cod{1:j-1}),2) + 1 : size(cat(2,cod{1:j}),2));
                        bempty{j} = zeros(p,size(cod{j},2));
                    end

                    % Calculates residuals of the full model
                    Xe = X - (cat(2,cod{:}) * Bmat');

                    % Calculates the residuals of the model without effect i
                    Xftmp = cod; Xftmp{i+1} = zeros(size(cod{i+1}));
                    Bf = B; Bf{i+1} = zeros(size(B{i+1}));
                    Ef = X - (cat(2,Xftmp{:}) * cat(2,Bf{:})');

                    % Calculates the sum of squares according to 
                    % Thiel et al. (2017)
                    ptest_res.ssq(k,i) = sum(Ef(:).^2) - sum(Xe(:).^2);

                else

                    Y = mfda_res.Y;

                    % Checks which coding is to be used to define the design
                    % matrices
                    if strcmp(Options.coding,'sumcod') == 1
                        [cod, codempty] = sumcoding(Y, Options); % sum coding of the input design matrix according to Thiel et al. (2017)
                    elseif strcmp(Options.coding,'wecod') == 1
                        [cod, codempty] = wecoding(Y, Options);  % weighted effect coding of the input design matrix according to te Grotenhuis et al. (2017)
                    else
                        error('The coding option defined does not exist or is not supported!');
                    end

                    % Random permutations of the coded effect
                    perms = randperm(n);
                    cod{i+1} = cod{i+1}(perms,:); % adds one to consider grand mean

                    % Calculates the model parameters in B
                    codcat = cat(2,cod{:});
                    Bmat = (pinv(codcat' * codcat) * codcat' * X)';

                    % Creates empty B and separates Bmat in array of cells (it is just useful for later calculations...)
                    bempty = {};
                    B = {};
                    for j = 1 : length(cod)
                        B{j} = Bmat(:,size(cat(2,cod{1:j-1}),2) + 1 : size(cat(2,cod{1:j}),2));
                        bempty{j} = zeros(p,size(cod{j},2));
                    end

                    % Calculates residuals of the full model
                    Xe = X - (cat(2,cod{:}) * Bmat');

                    % Calculates the residuals of the model without effect i
                    Xftmp = cod; Xftmp{i+1} = zeros(size(cod{i+1}));
                    Bf = B; Bf{i+1} = zeros(size(B{i+1}));
                    Ef = X - (cat(2,Xftmp{:}) * cat(2,Bf{:})');

                    % Calculates the sum of squares according to 
                    % Thiel et al. (2017)
                    ptest_res.ssq(k,i) = sum(Ef(:).^2) - sum(Xe(:).^2);


                end

            end

        end

    case 'rebalanced'

        for k = 1 : Options.nperms

            clc; fprintf('Performing permutations test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);

            % Loops over the effects
            for i = 1 : length(effects)

                if length(effects{i}) == 1 % main effects

                    Y = mfda_res.Y; % extracts true design matrix
                    Yo = Y(:,1:end ~= i); % design matrix of other factors than i
                    Youni = unique(Yo,'rows'); % unique combinations of other factors

                    % Creates randomized design matrix for factor i without
                    % affecting the other factors
                    for j = 1 : size(Youni)
                        idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                        perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                        Y(idx,i) = Y(idx(perms),i); % randomizes the design for factor i
                    end

                    % Coding of the design based on the permuted main effects
                    if strcmp(Options.coding,'sumcod') == 1

                        % Sum coding of the input design matrix according to
                        % Thiel et al. (2017)
                        [cod, codempty] = sumcoding(Y, Options);

                    elseif strcmp(Options.coding,'wecod') == 1

                        % Weighted-effect coding of the input design matrix according
                        % to Nieuwenhuis et al. (2017)
                        [cod, codempty] = wecoding(Y, Options);

                    end

                    % Coding of the rebalanced design based on the permuted 
                    % main effects
                    Options.decomp = 'glm';
                    [radecomp_res, Xbal, Ybal] = radecomp(X, Y, Options);
                    codbal = radecomp_res.cod;

                    % Calculates the model parameters in B (calibration)
                    codcat = cat(2,codbal{:});
                    Bmat = (pinv(codcat' * codcat) * codcat' * Xbal)';

                    % Creates empty B and separates Bmat in array of cells (it is just useful for later calculations...)
                    bempty = {};
                    B = {};
                    for j = 1 : length(cod)
                        B{j} = Bmat(:,size(cat(2,cod{1:j-1}),2) + 1 : size(cat(2,cod{1:j}),2));
                        bempty{j} = zeros(p,size(cod{j},2));
                    end

                    % Calculates residuals of the full model
                    Xe = X - (cat(2,cod{:}) * Bmat');

                    % Calculates the residuals of the model without effect i
                    % using the coding of the initial design and the B from
                    % the rebalanced design
                    Xftmp = cod; Xftmp{i+1} = zeros(size(cod{i+1}));
                    Bf = B; Bf{i+1} = zeros(size(B{i+1}));
                    Ef = X - (cat(2,Xftmp{:}) * cat(2,Bf{:})');

                    % Calculates the sum of squares according to 
                    % Thiel et al. (2017)
                    ptest_res.ssq(k,i) = sum(Ef(:).^2) - sum(Xe(:).^2);

                else

                    Y = mfda_res.Y;

                    % Checks which coding is to be used to define the design
                    % matrices
                    if strcmp(Options.coding,'sumcod') == 1
                        [cod, codempty] = sumcoding(Y, Options); % sum coding of the input design matrix according to Thiel et al. (2017)
                    elseif strcmp(Options.coding,'wecod') == 1
                        [cod, codempty] = wecoding(Y, Options);  % weighted effect coding of the input design matrix according to te Grotenhuis et al. (2017)
                    else
                        error('The coding option defined does not exist or is not supported!');
                    end

                    % Coding of the rebalanced design based on the permuted 
                    % main effects
                    Options.decomp = 'glm';
                    [radecomp_res, Xbal, Ybal] = radecomp(X, Y, Options);
                    codbal = radecomp_res.cod;
                    perms = randperm(size(Xbal,1));
                    codbal{i+1} = codbal{i+1}(perms,:); % adds one to consider grand mean

                    % Calculates the model parameters in B
                    codcat = cat(2,codbal{:});
                    Bmat = (pinv(codcat' * codcat) * codcat' * Xbal)';

                    % Creates empty B and separates Bmat in array of cells (it is just useful for later calculations...)
                    bempty = {};
                    B = {};
                    for j = 1 : length(cod)
                        B{j} = Bmat(:,size(cat(2,cod{1:j-1}),2) + 1 : size(cat(2,cod{1:j}),2));
                        bempty{j} = zeros(p,size(cod{j},2));
                    end

                    perms = randperm(n);
                    cod{i+1} = cod{i+1}(perms,:); % adds one to consider grand mean

                    % Calculates residuals of the full model
                    Xe = X - (cat(2,cod{:}) * Bmat');

                    % Calculates the residuals of the model without effect i
                    % using the coding of the initial design and the B from
                    % the rebalanced design
                    Xftmp = cod; Xftmp{i+1} = zeros(size(cod{i+1}));
                    Bf = B; Bf{i+1} = zeros(size(B{i+1}));
                    Ef = X - (cat(2,Xftmp{:}) * cat(2,Bf{:})');

                    % Calculates the sum of squares according to 
                    % Thiel et al. (2017)
                    ptest_res.ssq(k,i) = sum(Ef(:).^2) - sum(Xe(:).^2);


                end

            end

        end
end


% Retrieves the true/un-randomized ssq values
ptest_res.tssq = mfda_res.ssq; % true/un-randomized ssq values (cells 1 and 2 contain ssq of initial matrix X and grand mean)

% Total sum of squares for residuals are not kept
ptest_res.tssq = ptest_res.tssq(1:end-1);

% Calculates the p-values as the percentage of ssq based on
% permutations equal or greater than the true experimental ssq
ptest_res.pval = zeros(1,length(ptest_res.tssq));
for i = 1 : length(ptest_res.tssq)
    ptest_res.pval(1,i) = sum(ptest_res.ssq(:,i) >= ptest_res.tssq(1,i)) / Options.nperms;
end

end