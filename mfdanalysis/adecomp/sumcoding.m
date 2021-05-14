function [sumcod, sumcodempty] = sumcoding(Y, Options)
%
% The sumcoding function codes design matrices according to the
% General Linear Models methodology described in Thiel et al. (2017) using
% the sum coding or deviation coding strategy.
%
% See reference for more information.
% 
% Reference :
% ===========
% Thiel, M., Féraud, B., & Govaerts, B. (2017). ASCA+ and APCA+: 
% Extensions of ASCA and APCA in the analysis of unbalanced multifactorial 
% designs: Analyzing unbalanced multifactorial designs with ASCA+ and APCA+. 
% Journal of Chemometrics, 31(6), e2895. https://doi.org/10.1002/cem.2895
%
% Input arguments :
% =================
% Y : matrix of factors in the columns and levels for each sample 
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.interactions : the maximum number of interactions to consider
%
% 
% Output arguments :
% ==================
% sumcod : cell structure of the sum-coding design matrices
% sumcodempty : copy of sumcod where all values are set to zero
% 
% Usage :
% =======
% Options.interactions = 2;
% [sumcod] = sumcoding(Y, Options);
%
% Related function :
% ==================
% adecomp.m (performs ANOVA decomposition)
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

[n,q] = size(Y); % number of experiments x factors

sumcod = {}; % allocates space for the sum coding
sumcod{1,1} = ones(n,1); % sum coding for the grand mean matrix

% Main effects
for i = 1 : q
    sumcod{1,i+1} = yformatconv(Y(:,i),'binarymat'); % converts Y integer vector to binary matrix
    sumcod{1,i+1}(sumcod{1,i+1}(:,end) == 1,:) = -1; % identifies last level position in rows and replaces the whole rows with -1
    sumcod{1,i+1}(:,end) = []; % deletes the last level column
end

% Sum coding for interactions only if Options.interactions >= 2
if Options.interactions >= 2
    
    % Defines interactions matrix according to Options.interactions
    % combs contains all possible interactions according to Options.interactions
    % Each cell in combs contains a combination of the main factors indices
    combs = {};
    for i = 2 : Options.interactions
        combs = [combs; num2cell(nchoosek(1:q, i),2)];
    end
    
    % Loops over the interactions in combs
    for j = 1 : length(combs)
        
        pos = q + j + 1; % defines position taking into account main effects and grand mean design matrices
        
        % Loops used to calculate pairwise Kronecker products between matrices
        for k = 1 : (length(unique(combs{j})) - 1) 
            
            if k == 1 % the first Kronecker product between two matrices
                sumcod{1,pos} = kron(sumcod{1,combs{j}(1)+1},sumcod{1,combs{j}(2)+1});
                sumcod{1,pos} = sumcod{1,pos}(1:n+1:end,:);
            else
                sumcod{1,pos} = kron(sumcod{1,pos},sumcod{1,combs{j}(k+1)+1});
                sumcod{1,pos} = sumcod{1,pos}(1:n+1:end,:);
            end
        
        end
    end

end

% sumcodeempty as an empty cell structure is useful in the main functions
sumcodempty = {};
for i = 1 : length(sumcod)
    sumcodempty{i} = zeros(size(sumcod{i}));
end

end