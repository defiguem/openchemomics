function [ptest_res] = mfdptest(mfda_res)
%
% The mfdptest function performs permutation tests in order to
% determine main effects/interactions significance. To do so, repeated
% randomizations of the rows in the design matrix are performed before
% calculating matrices averages.
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
%
% =========================================================================

%% Fail-safe section

[n,p] = size(mfda_res.avgmat{1}); % size of X
[~,q] = size(mfda_res.Y); % size of Y

% Retrieves field from mfda_res structure
Options = mfda_res.Options;

% Checks if the number of permutation tests to perform by effect was chosen
if isfield(Options,'nperms') == 0
    Options.nperms = 1000;
end

if strcmp(Options.decomp,'glm') == 1 && isfield(Options,'coding') == 0
    Options.coding = 'sumcod';
end

%% Permutation tests

% Note that in mfda_res, the ANOVA matrix after deflation of the grand mean
% matrix is stored in : mfda_res.anovamat{1,2};

ptest_res.ssq = zeros(Options.nperms, length(mfda_res.matid)-2); % allocates space for ssq of permuted levels

% Following line makes only a difference if RASCA was used
mfda_res.Y = cat(2,mfda_res.Ys{1:q});

switch Options.decomp
    
    case {'classical','res_classical'}
        
        for k = 1 : Options.nperms
            
            % clc; fprintf('Performing permutation test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);
            
            % Main effects
            for i = 1 : q
                
                % Creates randomized design matrix for factor i
                Y = mfda_res.Y; % extracts true design matrix
                Yo = Y(:,1:end~=i); % Design matrix of other factors than i
                Youni = unique(Yo,'rows'); % unique combinations of other factors
                
                for j = 1 : size(Youni)
                    idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                    perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                    Y(idx,i) = Y(idx(perms),i); % randomizes the design for factor i
                end
                
                nl = max(Y(:,i)); % number of levels for factor Y(:,i)
                
                % Creates matrix of averages for each level j of factor i
                avgmat = zeros(n,p); % initializes average matrix
                for j = 1 : nl
                    idx = find(Y(:,i) == j); % indices of samples at factor i and level j
                    %levmean = mean(mfda_res.anovamat{1,i+1}(idx,:),1); % mean of factor i at level j
                    levmean = mean(mfda_res.anovamat{1,2}(idx,:),1); % mean of factor i at level j
                    avgmat(idx,:) = repmat(levmean, length(idx),1);
                end
                
                % Calculates the sum of squares for permutated level assignments
                ptest_res.ssq(k,i) = sum(avgmat(:).^2);
                
            end
            
            % Interactions
            if Options.interactions >= 2
                
                % Defines interactions matrix according to Options.interactions
                combs = {};
                for i = 2 : Options.interactions
                    combs = [combs; num2cell(nchoosek(1:q, i),2)];
                end
                
                tab = q + 1; % adds one position to last table
                
                % For interactions, the whole matrix is randomized
                perms = randperm(n); % random permutations
                Y = mfda_res.Y(perms,:); % randomizes levels assignment
                
                for i = 1 : length(combs)
                    
                    Ycomb = cellstr(num2str(Y(:,combs{i}))); % interaction of factors
                    uni = unique(Ycomb); % unique combinations of levels in factors
                    
                    % Creates matrix of averages for each level j of factor i
                    avgmat = zeros(n,p); % initializes average matrix
                    for j = 1 : length(uni)
                        idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of samples at factor i and level j
                        %levmean = mean(mfda_res.anovamat{1,tab+1}(idx,:),1); % mean of factor i at level j
                        levmean = mean(mfda_res.anovamat{1,2+size(Y,2)}(idx,:),1); % mean of factor i at level j
                        avgmat(idx,:) = repmat(levmean, length(idx),1);
                    end
                    
                    % Calculates the sum of squares for permutated level assignments
                    ptest_res.ssq(k,tab) = sum(avgmat(:).^2); % sum of squares
                    
                    tab = tab + 1; % updates table position
                    
                end
                
            end
            
        end
        
    case 'glm'
        
        for k = 1 : Options.nperms
            
            %clc; fprintf('Performing permutation test, please wait...\n Iteration : %d / %d \n', k, Options.nperms);
            
            % Calculation of the main factor effects : iterates over columns of Y
            for i = 1 : q
                
                % Creates randomized design matrix for factor i
                Y = mfda_res.Y; % extracts true design matrix
                Yo = Y(:,1:end~=i); % Design matrix of other factors than i
                Youni = unique(Yo,'rows'); % unique combinations of other factors
                
                for j = 1 : size(Youni)
                    idx = find( ismember(Yo,Youni(j,:),'rows') == 1 ); % finds position of each unique combination of other factors
                    perms  = randperm(length(idx)); % creates randomized indices for the length of idx
                    Y(idx,i) = Y(idx(perms),i); % randomizes the design for factor i
                end
                
                if strcmp(Options.coding,'sumcod') == 1
                    % Sum coding of the input design matrix according to Thiel et al. (2017)
                    [cod, codempty] = sumcoding(Y, Options);
                elseif strcmp(Options.coding,'wecod') == 1
                    [cod, codempty] = wecoding(Y, Options);
                else
                    error('The coding option defined does not exist or is not supported!');
                end
                
                % Creates empty B (it is just useful for later calculations...)
                bempty = {};
                for j = 1 : length(cod)
                    bempty{j} = zeros(p,size(cod{j},2));
                end
                
                B = bempty; B{i+1} = (pinv(cod{i+1}' * cod{i+1}) * cod{i+1}' * mfda_res.anovamat{1,i+1})';
                Xf = codempty; Xf{i+1} = cod{i+1};
                
                avgmat = cat(2,Xf{:}) * cat(2,B{:})';
                
                % Calculates the sum of squares for permutated level assignments
                ptest_res.ssq(k,i) = sum(avgmat(:).^2);
                
            end
            
            % Calculates interactions only if Options.interactions >= 2
            if Options.interactions >= 2
                
                % For interactions, the whole matrix is randomized
                perms = randperm(n); % random permutations
                Y = mfda_res.Y(perms,:); % randomizes levels assignment
                
                if strcmp(Options.coding,'sumcod') == 1
                    % Sum coding of the input design matrix according to Thiel et al. (2017)
                    [cod, codempty] = sumcoding(Y, Options);
                elseif strcmp(Options.coding,'wecod') == 1
                    [cod, codempty] = wecoding(Y, Options);
                end
                
                % Defines interactions matrix according to Options.interactions
                combs = {};
                for i = 2 : Options.interactions
                    combs = [combs; num2cell(nchoosek(1:q, i),2)];
                end
                
                tab = q + 1;
                
                for i = 1 : length(combs)
                    
                    B = bempty; B{tab} = (pinv(cod{tab}' * cod{tab}) * cod{tab}' * mfda_res.anovamat{1,tab+1})';
                    Xf = codempty; Xf{tab} = cod{tab};
                    
                    avgmat = cat(2,Xf{:}) * cat(2,B{:})';
                    
                    % Calculates the sum of squares for permutated level assignments
                    ptest_res.ssq(k,tab) = sum(avgmat(:).^2); % sum of squares
                    
                    tab = tab + 1; % updates table position
                    
                end
                
            end
            
        end
end

ptest_res.tssq = cell2mat(mfda_res.ssq(3:end)); % true/un-randomized ssq values (cells 1 and 2 contain ssq of initial matrix X and grand mean)
ptest_res.tssq = ptest_res.tssq(1:end-1); % total sum of squares for residuals are not kept

ptest_res.pval = zeros(1,length(ptest_res.tssq));
for i = 1 : length(ptest_res.tssq)
    ptest_res.pval(1,i) = sum(ptest_res.ssq(:,i) >= ptest_res.tssq(1,i)) / Options.nperms;
end

end