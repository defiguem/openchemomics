function [ptest_res] = mfdptest(mfda_res)
%
% The mfdptest function performs permutation tests in order to
% determine main effects/interactions significance. To do so, repeated
% randomizations of the rows in the design matrix are performed before
% calculating effect matrices.
%
% Whereas for interaction the rows of the design matrix are permutated
% without restrictions, for the main effects there is one. Permutations
% for one factor are only allowed within the levels of other factors. For
% more details on permutation tests, see Anderson et Braak, 2003.
%
% Moreover, permutations for the interaction are made on the matrix of
% residuals after removing all main effects.
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
% Input arguments :
% =================
% mfda_res : output argument of the adecomp.m function or of any other
% supported function performing analysis of multifactorial designs such as
% apca.m, asca.m, acomdim.m, etc.
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
% Related functions :
% ===================
% adecomp.m (performs multivariate ANOVA decompostion)
% apca.m
% asca.m
% acomdim.m
% 
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
% 2021.08.05 : adapted permutations tests to multifactorial design analysis
%              based on a resampling strategy
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

% Retrieves the successively decomposed ANOVA matrices depending on the use
% of a resampling strategy or not
if strcmp(Options.decomp,'classical') == 1 || strcmp(Options.decomp,'glm') == 1 
    anovamat = mfda_res.anovamat;
elseif strcmp(Options.decomp,'res_classical') == 1 || strcmp(Options.decomp,'res_glm') == 1 % based on a resampling strategy
    anovamat = mfda_res.panovamat;
end

[n,p] = size(anovamat{1}); % size of X
[~,q] = size(mfda_res.Y); % size of Y

% Allocates space for ssq of permuted levels
ptest_res.ssq = zeros(Options.nperms, length(mfda_res.matid)-2); 

switch Options.decomp
    
    case {'classical','res_classical'}
        
        for k = 1 : Options.nperms
            
            % clc; fprintf('Performing permutation test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);
            
            % Main effects
            for i = 1 : q
                
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
                
                nl = max(Y(:,i)); % number of levels for factor Y(:,i)
                
                % Creates matrix of averages for each level j of factor i
                avgmat = zeros(n,p); % initializes average matrix
                for j = 1 : nl
                    idx = find(Y(:,i) == j); % indices of observations of factor i at level j
                    levmean = mean(anovamat{1,1}(idx,:),1); % mean of observations of factor i at level j
                    avgmat(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                end
                
                % Calculates the sum of squares for permutated level assignments
                ptest_res.ssq(k,i) = sum(avgmat(:).^2);
                
            end
            
            % Interactions
            if Options.interactions >= 2
                
                % Defines interactions matrix according to Options.interactions
                % Each cell contains a possible combination of factors
                combs = {};
                for i = 2 : Options.interactions
                    combs = [combs; num2cell(nchoosek(1:q, i),2)];
                end
                
                % Adds one position to last table to take into account the
                % main effects
                tab = q + 1;
                
                % For interactions, the whole matrix is randomized
                perms = randperm(n); % random permutations
                Y = mfda_res.Y(perms,:); % randomizes levels assignment
                
                % Loops over the possible combinations of levels (can be
                % interactions with any number of factors)
                for i = 1 : length(combs)
                    
                    Ycomb = cellstr(num2str(Y(:,combs{i}))); % interaction of factors in combs{i}
                    uni = unique(Ycomb); % defines unique combinations of levels in factors combs{i}
                    
                    % Creates matrix of averages for each unique combination 
                    % of factors combs{i} at levels uni{j}
                    avgmat = zeros(n,p); % initializes average matrix
                    for j = 1 : length(uni)
                        idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of observations of factors combs{i} at levels uni{j}
                        levmean = mean(anovamat{1,2+size(Y,2)}(idx,:),1); % mean of observations of factors combs{i} at levels uni{j}
                        avgmat(idx,:) = repmat(levmean, length(idx),1); % replicates the mean vector in the effect matrix
                    end
                    
                    % Calculates the sum of squares for permutated level assignments
                    ptest_res.ssq(k,tab) = sum(avgmat(:).^2); % sum of squares
                    
                    tab = tab + 1; % updates table position
                    
                end
                
            end
            
        end
        
    case {'glm','res_glm'}
        
        for k = 1 : Options.nperms
            
            % clc; fprintf('Performing permutation test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);
            
            % Calculation of the main factor effects : iterates over columns of Y
            for i = 1 : q
                
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
                
                % Checks which coding is to be used to define the design
                % matrices
                if strcmp(Options.coding,'sumcod') == 1
                    [cod, codempty] = sumcoding(Y, Options); % sum coding of the input design matrix according to Thiel et al. (2017)
                elseif strcmp(Options.coding,'wecod') == 1
                    [cod, codempty] = wecoding(Y, Options);  % weighted effect coding of the input design matrix according to te Grotenhuis et al. (2017)
                else
                    error('The coding option defined does not exist or is not supported!');
                end
                
                % Creates empty B (it is just useful for later calculations...)
                bempty = {};
                for j = 1 : length(cod)
                    bempty{j} = zeros(p,size(cod{j},2));
                end
                
                % Copies the empty B array of cells and calculates the 
                % parameters matrix for factor i
                B = bempty;
                B{i+1} = (pinv(cod{i+1}' * cod{i+1}) * cod{i+1}' * anovamat{1,1})'; 
                
                % Copies the empty array of cells for the design matrices
                % and adds the design matrix of factor i in Xf, the rest 
                % are zeros
                Xf = codempty;
                Xf{i+1} = cod{i+1}; 
                
                % Calculates the effect matrix of factor i 
                avgmat = cat(2,Xf{:}) * cat(2,B{:})'; 
                
                % Calculates the sum of squares for permutated level assignments
                ptest_res.ssq(k,i) = sum(avgmat(:).^2);
                
            end
            
            % Calculates interactions only if Options.interactions >= 2
            if Options.interactions >= 2
                
                % For interactions, the whole matrix is randomized
                perms = randperm(n); % random permutations
                Y = mfda_res.Y(perms,:); % randomizes levels assignment
                
                % Checks which coding is to be used to define the design
                % matrices
                if strcmp(Options.coding,'sumcod') == 1
                    [cod, codempty] = sumcoding(Y, Options); % sum coding of the input design matrix according to Thiel et al. (2017)
                elseif strcmp(Options.coding,'wecod') == 1
                    [cod, codempty] = wecoding(Y, Options);  % weighted effect coding of the input design matrix according to te Grotenhuis et al. (2017)
                else
                    error('The coding option defined does not exist or is not supported!');
                end
                
                % Defines interactions matrix according to Options.interactions
                % Each cell contains a possible combination of factors
                combs = {};
                for i = 2 : Options.interactions
                    combs = [combs; num2cell(nchoosek(1:q, i),2)];
                end
                
                % Adds one position to last table to take into account the
                % main effects
                tab = q + 1;
                
                % Loops over the possible combinations of levels (can be
                % interactions with any number of factors)
                for i = 1 : length(combs)
                    
                    % Copies the empty B array of cells and calculates the 
                    % parameters matrix
                    B = bempty; 
                    B{tab} = (pinv(cod{tab}' * cod{tab}) * cod{tab}' * anovamat{1,1})';
                    
                    % Copies the empty array of cells for the design matrices
                    % and adds the design matrix of the interaction in 
                    % combs{i} in Xf, the rest are zeros
                    Xf = codempty;
                    Xf{tab} = cod{tab};
                    
                    % Calculates the effect matrix of the interaction in combs{i}
                    avgmat = cat(2,Xf{:}) * cat(2,B{:})';
                    
                    % Calculates the sum of squares for permutated level assignments
                    ptest_res.ssq(k,tab) = sum(avgmat(:).^2); % sum of squares
                    
                    tab = tab + 1; % updates table position
                    
                end
                
            end
            
        end
end

% Retrieves the true/un-randomized ssq values depending on the use
% of a resampling strategy or not
if strcmp(Options.decomp,'classical') == 1 || strcmp(Options.decomp,'glm') == 1
    ptest_res.tssq = cell2mat(mfda_res.ssq(3:end)); % true/un-randomized ssq values (cells 1 and 2 contain ssq of initial matrix X and grand mean)
elseif strcmp(Options.decomp,'res_classical') == 1 || strcmp(Options.decomp,'res_glm') == 1
    ptest_res.tssq = cell2mat(mfda_res.pssq(3:end)); % true/un-randomized ssq values (cells 1 and 2 contain ssq of initial matrix X and grand mean)
end

% Total sum of squares for residuals are not kept
ptest_res.tssq = ptest_res.tssq(1:end-1); 

% Calculates the p-values as the percentage of ssq based on
% permutations equal or greater than the true experimental ssq
ptest_res.pval = zeros(1,length(ptest_res.tssq));
for i = 1 : length(ptest_res.tssq)
    ptest_res.pval(1,i) = sum(ptest_res.ssq(:,i) >= ptest_res.tssq(1,i)) / Options.nperms;
end

end