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
% multivariate, are investigated to assert whether such variables have an
% effect on the system under study.
%
% The decomposition can be done in two ways :
% (1) classical multivariate ANOVA decomposition
% into matrix effects for factors/interactions by
% successively extracting the grand mean matrix, main factor effects and
% interactions based on the means of factors, levels, and interactions;
% (2) multivariate ANOVA decomposition based on
% the General Linear Models (GLM) presented by Thiel et al. (2017).
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
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.decomp : ANOVA decomposition method 'classical', 'glm',
%       'res_classical'
%   Options.coding : if Options.decomp is 'glm', the coding method of the
%       design matrices can be chosen as 'sumcod' or 'wecod'
%   Options.interactions : the maximum number of interactions to consider
%
% Output arguments :
% ==================
% adecomp_res : structure containing results of the ANOVA decomposition
%   adecomp_res.anovamat : ANOVA matrices contained in cells (the last cell
%       contains the residuals)
%   adecomp_res.avgmat : cell structure of factor/interation effects
%       (see adecomp_res.matid for information on factor/interactions)
%   adecomp_res.matid : identity of ANOVA matrices (grand mean, effects,
%       interactions and residual pure error)
%   adecomp_res.ssq : sum of squares for each effect based on adecomp_res.avgmat
%   adecomp_res.ssqvarexp : variance explained by the sum of squares
%   adecomp_res.Ys : stores the design matrices of factors and interactions
%   adecomp_res.Xs : stores the effect matrix to which the pure error 
%       was added
%   adecomp_res.Ysm : mean effect vector for each factor and interaction
%
%   adecomp_res.Options : Options used to perform ANOVA decomposition
%
% Note that for the explained variance, the reference total sum of squares
% is calculated based on the sum of squares of the first deflated matrix in
% position adecomp_res.anovamat{1,2} corresponding to the mean centered
% form of the input argument X, in other words after grand mean
% substraction.
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
    Options.interactions = 2;
elseif isfield(Options,'interactions') == 1 && Options.interactions > q
    Options.interactions = q; % maximum number of interactions cannot be > q
end

% Checks if the ANOVA decomposition method was specified
if isfield(Options,'decomp') == 0
    Options.decomp = 'classical';
end

% Checks if the ANOVA decomposition method was specified
if isfield(Options,'decomp') == 1 && strcmp(Options.decomp,'glm') == 1 && isfield(Options,'coding') == 0
    Options.coding = 'sumcod';
end

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : q
    Y(:,i) = yformatconv(Y(:,i),'intvec');
end

%% Multivariate ANOVA decomposition

% The input argument X is placed as the first (not deflated) cell
adecomp_res.anovamat{1,1} = X;

% Switches cases depending on the multivariate ANOVA decomposition method
switch Options.decomp
    
    case {'classical','res_classical'}
        
        % Grand mean matrix : initialization of decomposition into ANOVA matrices
        adecomp_res.avgmat{1,1} = repmat(mean(X,1), n,1); % grand mean matrix
        adecomp_res.matid{1,1} = 'Mean'; % grand mean matrix as index 0
        
        % Deflates input argument X by substracting grand mean matrix
        adecomp_res.anovamat{1,2} = X - adecomp_res.avgmat{1,1};
        
        % Calculates sum of squares of the initial matrix
        adecomp_res.ssq{1,1} = sum(adecomp_res.anovamat{1,1}(:).^2); % total sum of squares
        
        % Calculates sum of squares of the grand mean matrix
        adecomp_res.ssq{1,2} = sum(adecomp_res.avgmat{1,1}(:).^2); % total sum of squares
        
        % Calculation of the main factor effects : iterates over columns of Y
        for i = 1 : q
            
            % Stores factor/level classes
            adecomp_res.Ys{1,i} = Y(:,i);
            
            % Stores levels for factor Y(:,i)
            levels = unique(Y(:,i)); 
            
            % Creates matrix of averages for each level j of factor i
            adecomp_res.anovamat{1,i+2} = zeros(n,p); % allocates space
            adecomp_res.avgmat{1,i+1} = zeros(n,p); % allocates space
            
            % Loops over the levels of a given factor
            for j = 1 : length(levels)
                idx = find(Y(:,i) == levels(j)); % indices of samples at factor i and level j
                levmean = mean(adecomp_res.anovamat{1,i+1}(idx,:),1); % mean of factor i at level j
                adecomp_res.avgmat{1,i+1}(idx,:) = repmat(levmean, length(idx),1); % places level means
            end
            
            % Replaces inplace the matrix of averages by decomposition from previous one
            adecomp_res.anovamat{1,i+2} = adecomp_res.anovamat{1,i+1} - adecomp_res.avgmat{1,i+1};
            adecomp_res.matid{1,i+1} = i; % factor index associated to Y column
            
            % Calculates sum of squares of average matrix for factor i
            adecomp_res.ssq{1,i+2} = sum(adecomp_res.avgmat{1,i+1}(:).^2); % sum of squares
            
            % Calculates percentage of explained variance of average matrix for factor i
            adecomp_res.ssqvarexp{1,i} = (adecomp_res.ssq{1,i+2} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2})) * 100; % explained variance
            
        end
        
        % Calculates interactions only if Options.interactions >= 2
        if Options.interactions >= 2
            
            % Defines interactions matrix according to Options.interactions
            combs = {};
            for i = 2 : Options.interactions
                combs = [combs; num2cell(nchoosek(1:q, i),2)];
            end
            
            tab = length(adecomp_res.avgmat) + 1; % Need to add one position to last table
            
            % Loops over combinations of factors (interactions)
            for i = 1 : length(combs)
                
                % Combination of factors (interactions)
                adecomp_res.Ys{1,i+q} = Y(:,combs{i});
                Ycomb = cellstr(num2str(adecomp_res.Ys{1,i+q}));
                
                % Unique combinations of levels in factors
                uni = unique(Ycomb); 
                
                % Creates matrix of averages for each unique combinations 
                % of levels in factors
                adecomp_res.anovamat{1,tab+1} = zeros(n,p);
                adecomp_res.avgmat{1,tab} = zeros(n,p);
                
                % Loops over combinations of levels in factors
                for j = 1 : length(uni)
                    idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of samples at factor i and level j
                    levmean = mean(adecomp_res.anovamat{1,tab}(idx,:),1); % mean of factor i at level j
                    adecomp_res.avgmat{1,tab}(idx,:) = repmat(levmean, length(idx),1);
                end
                
                % Replaces in place the matrix of averages by decomposition from previous one
                adecomp_res.anovamat{1,tab+1} = adecomp_res.anovamat{1,tab} - adecomp_res.avgmat{1,tab};
                adecomp_res.matid{1,tab} = combs{i}; % factor index associated to Y column
                
                % Calculates sum of squares of average matrix for combinations of levels in factors' interactions
                adecomp_res.ssq{1,tab+1} = sum(adecomp_res.avgmat{1,tab}(:).^2); % sum of squares
                
                % Calculates percentage of explained variance for combinations of levels in factors' interactions
                adecomp_res.ssqvarexp{1,tab-1} = (adecomp_res.ssq{1,tab+1} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2})) * 100; % explained variance
                
                tab = tab + 1; % updates table position
                
            end
            
        end

    case {'glm','res_glm'}
        
        % Either codes design matrices using the sum coding or the
        % weighted-effect coding
        if strcmp(Options.coding,'sumcod') == 1
            
            % Sum coding of the input design matrix according to 
            % Thiel et al. (2017)
            [cod, codempty] = sumcoding(Y, Options);
            adecomp_res.sumcod = cod;
            
        elseif strcmp(Options.coding,'wecod') == 1
            
            % Weighted-effect coding of the input design matrix according
            % to Nieuwenhuis et al. (2017)
            [cod, codempty] = wecoding(Y, Options);
            adecomp_res.wecod = cod;
            
        end
        
        % Creates empty B (it is just useful for later calculations...)
        bempty = {};
        for i = 1 : length(cod)
            bempty{i} = zeros(p,size(cod{i},2));
        end
        
        % Grand mean matrix : initialization of decomposition into GLM matrices
        B = bempty;
        B{1} = (pinv(cod{1}' * cod{1}) * cod{1}' * X)';
        Xf = codempty;
        Xf{1} = cod{1};
        
        adecomp_res.avgmat{1,1} = cat(2,Xf{:}) * cat(2,B{:})'; % grand mean matrix 
        
        adecomp_res.anovamat{1,2} = X - adecomp_res.avgmat{1,1}; % centered X
        adecomp_res.matid{1,1} = 'Mean'; % grand mean matrix as index 0
        
        % Calculates sum of squares of the initial matrix
        adecomp_res.ssq{1,1} = sum(adecomp_res.anovamat{1,1}(:).^2); % total sum of squares
        
        % Calculates sum of squares of the grand mean matrix
        adecomp_res.ssq{1,2} = sum(adecomp_res.avgmat{1,1}(:).^2); % total sum of squares
        
        % Calculation of the main factor effects : iterates over columns of Y
        for i = 1 : q
            
            % Stores factor/level classes
            adecomp_res.Ys{1,i} = Y(:,i); 

            B = bempty;
            B{i+1} = (pinv(cod{i+1}' * cod{i+1}) * cod{i+1}' * adecomp_res.anovamat{1,i+1})';
            
            Xf = codempty;
            Xf{i+1} = cod{i+1};
            
            adecomp_res.avgmat{1,i+1} = cat(2,Xf{:}) * cat(2,B{:})';
            
            % Replaces inplace the matrix of averages by decomposition from previous one
            adecomp_res.anovamat{1,i+2} = adecomp_res.anovamat{1,i+1} - adecomp_res.avgmat{1,i+1};
            adecomp_res.matid{1,i+1} = i; % factor index associated to Y column
            
            % Calculates sum of squares of average matrix for factor i
            adecomp_res.ssq{1,i+2} = sum(adecomp_res.avgmat{1,i+1}(:).^2); % sum of squares
            
            % Calculates percentage of explained variance of average matrix for factor i
            adecomp_res.ssqvarexp{1,i} = (adecomp_res.ssq{1,i+2} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2})) * 100; % explained variance
            
        end
        
        % Calculates interactions only if Options.interactions >= 2
        if Options.interactions >= 2
            
            % Defines interactions matrix according to Options.interactions
            combs = {};
            for i = 2 : Options.interactions
                combs = [combs; num2cell(nchoosek(1:q, i),2)];
            end
            
            tab = length(adecomp_res.avgmat)+1; % Need to add one position to last table
            
            for i = 1 : length(combs)
                
                % Combination of factors (interactions)
                adecomp_res.Ys{1,i+q} = Y(:,combs{i});
                Ycomb = cellstr(num2str(adecomp_res.Ys{1,i+q}));
                
                B = bempty;
                B{tab} = (pinv(cod{tab}' * cod{tab}) * cod{tab}' * adecomp_res.anovamat{1,tab})';
                
                Xf = codempty;
                Xf{tab} = cod{tab};
                
                adecomp_res.avgmat{1,tab} = cat(2,Xf{:}) * cat(2,B{:})';
                
                % Replaces in place the matrix of averages by decomposition from previous one
                adecomp_res.anovamat{1,tab+1} = adecomp_res.anovamat{1,tab} - adecomp_res.avgmat{1,tab};
                adecomp_res.matid{1,tab} = combs{i}; % factor index associated to Y column
                
                % Calculates sum of squares of average matrix for combinations of levels in factors' interactions
                adecomp_res.ssq{1,tab+1} = sum(adecomp_res.avgmat{1,tab}(:).^2); % sum of squares
                
                % Calculates percentage of explained variance for combinations of levels in factors' interactions
                adecomp_res.ssqvarexp{1,tab-1} = (adecomp_res.ssq{1,tab+1} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2})) * 100; % explained variance
                
                tab = tab + 1; % updates table position
                
            end
        end
end

% Number of effects (grand mean, factors, interactions)
ntab = size(adecomp_res.avgmat,2);

% Loops over the effects, except for the grand mean
Xs = cell(1, ntab);
for i = 1 : ntab - 1
    Xs{1,i} = adecomp_res.avgmat{1,i+1} + adecomp_res.anovamat{1,end};
end
  
% Adds the residuals matrix as a last table and additional information 
Xs{1,ntab} = adecomp_res.anovamat{1,end};
adecomp_res.matid{1,ntab+1} = 'Residuals'; 
adecomp_res.ssq{1,ntab+2} = sum(Xs{1,ntab}(:).^2); 
adecomp_res.ssqvarexp{1,ntab} = (adecomp_res.ssq{1,ntab+2} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2})) * 100; 

% Calculation of the type III sum of squares
for i = 1 : ntab
    
    if i ~= ntab
        Ef = (adecomp_res.anovamat{end} + adecomp_res.avgmat{i+1}).^2;
        adecomp_res.ssqvarexp{2,i} = (sum(Ef(:)) - adecomp_res.ssq{1,end}) / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2}) * 100;
    else
        adecomp_res.ssqvarexp{2,i} = adecomp_res.ssq{1,end} / (adecomp_res.ssq{1,1}-adecomp_res.ssq{1,2}) * 100;
    end
    
end

% Place effect mean matrices + residuals in output
adecomp_res.Xs = Xs;

% Defines what is the mean vector for each factor and interaction
adecomp_res.Ysm = cell(1,length(adecomp_res.Ys));
for i = 1 : length(adecomp_res.Ys)
    
    adecomp_res.Ysm{1,i}{1} = unique(adecomp_res.Ys{i},'rows');
    
    for j = 1 : size(adecomp_res.Ysm{1,i}{1},1)
        
        [~,idx] = ismember(adecomp_res.Ys{i},adecomp_res.Ysm{1,i}{1}(j,:), 'rows');
        adecomp_res.Ysm{1,i}{2}(j,:) = mean(adecomp_res.avgmat{i+1}(idx == 1, :),1);
        
    end

end

end
