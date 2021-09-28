function [wecodfull, wecodempty] = wecoding(Y, Options)
%
% The wecoding function defines coding matrices according to the
% General Linear Models methodology with the specific goal of tackling the
% unbalancedness of designs of experiments. "wecoding" stands for "weighted 
% effect coding". Instead of considering deviations of main effects and 
% interactions from the global mean, wecoding uses weighted means to 
% account for prevalences in the levels of the different factors.
% 
% References :
% ============
% Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., Schmidt-Catran, 
% A., & Konig, R. (2017). When size matters : Advantages of weighted effect 
% coding in observational studies. International Journal of Public Health, 
% 62(1), 163-167. https://doi.org/10.1007/s00038-016-0901-1
%
% te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., 
% Schmidt-Catran, A., & Konig, R. (2017). A novel method for modelling 
% interaction between categorical variables. International Journal of 
% Public Health, 62(3), 427-431. https://doi.org/10.1007/s00038-016-0902-0
%
% URL : http://www.ru.nl/sociology/mt/wec/downloads/
%   Check the "Explaining the weighted coded interactions in Excel" section
%   to download Excel files with weighted effect coding rules for main
%   effects and interactions.
%
% Input arguments :
% =================
% Y : matrix of factors in the columns and levels for each sample 
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.interactions : the maximum number of interactions to consider
% 
% Output arguments :
% ==================
% wecod : cell structure of the wecoding results
%
% wecodempty : copy of wecod where all values are set to zero
% 
% Usage :
% =======
% Options.interactions = 2;
% [wecod] = wecoding(Y, Options);
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

% Finds the unique combinations of levels between the factors
Ycomb = unique(Y,'rows','legacy');

% Finds the maximum number of levels per factor
maxfaclvl = zeros(1,q);
keptlvls = cell(1,q);
for i = 1 : q
    maxfaclvl(i) = max(Ycomb(:,i));
end

% Defines by default the levels to omit as the last levels of each factor
% Add option in case the omitted levels should be different
omit = maxfaclvl;

% Finds the levels to be kept per factor
for i = 1 : q
    keptlvls{i} = 1 : maxfaclvl(i);
    keptlvls{i} = keptlvls{i}(1:end~=omit(i));
end

% Calculates the frequencies of each combination of levels in a table
freqtab = zeros(maxfaclvl);
for i = 1 : size(Ycomb,1)
    idx = num2cell(Ycomb(i,:));
    freqtab(idx{:}) = sum(ismember(Y,Ycomb(i,:),'rows'));
end

% Allocates space for the WE coding
wecod = {}; 
wecod{1,1} = ones(size(Ycomb,1),1); % we coding for the grand mean matrix
    
% Design matrices of main effects 
for i = 1 : q
    wecod{1,i+1} = yformatconv(Ycomb(:,i),'binarymat'); % converts Y integer vector to binary matrix
    freq = histcounts(Y(:,i)); % calculate the frequencies of each level in factor i    
    wecod{1,i+1}(wecod{1,i+1}(:,omit(i)) == 1,1:end~=omit(i)) = repmat(-1 * freq(1:end~=omit(i))./freq(omit(i)),[sum(wecod{1,i+1}(:,omit(i)) == 1),1]); 
    wecod{1,i+1} = wecod{1,i+1}(:,1:end~=omit(i)); % deletes the last level column
end

% Weighted effect coding for interactions only if Options.interactions == 2
if Options.interactions == 2
    
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
        
        wecod{1,pos} = kron(wecod{1,combs{j}(1)+1},wecod{1,combs{j}(2)+1});
        wecod{1,pos} = wecod{1,pos}(1:size(Ycomb,1)+1:end,:);
        [v1,v2] = ndgrid(1:size(wecod{1,combs{j}(1)+1}',1),1:size(wecod{1,combs{j}(2)+1},2));
        lvlscomb = sortrows([reshape(v1,[],1),reshape(v2,[],1)],1);
        
        for k = 1 : size(lvlscomb,1)
            
            prodsigns = [sign(wecod{1,combs{j}(1)+1}(:,lvlscomb(k,1))),sign(wecod{1,combs{j}(2)+1}(:,lvlscomb(k,2)))];
            
            % Condition 1 : whenever the product of the two weighted effect 
            % coded dummies are >=0 and none of the individual effects has
            % a negative values, the coding on the interaction is simply 
            % their product.
            idx = ( (prodsigns(:,1) .* prodsigns(:,2)) >= 0 & sum(prodsigns,2) ~= -2);
            wecod{1,pos}(idx == 1,k) = wecod{1,combs{j}(1)+1}(idx,lvlscomb(k,1)) .* wecod{1,combs{j}(2)+1}(idx,lvlscomb(k,2));
            
            % Condition 2 : it is an exception to condition 1. If there are
            % negative values for both weighted effect coded dummies, the 
            % coding is a positive fraction where the numerator is the
            % number of observations in the cell for which the two product 
            % terms are indicators divided by the number of observations in
            % the omitted category
            idx = sum(prodsigns,2) == -2;
            tabidx = num2cell([keptlvls{combs{j}(1)}(lvlscomb(k,1)),keptlvls{combs{j}(2)}(lvlscomb(k,2))]);
            omitidx = num2cell(omit);
            wecod{1,pos}(idx == 1,k) = freqtab(tabidx{:}) ./ freqtab(omitidx{:});
            
            % Condition 3 : whenever the product of two weighted  effect 
            % coded dummies are < 0, a negative fraction is inserted as
            % the number of observations in the cell for which the two 
            % product terms are indicators divided by the number of 
            % observations in the cell to which the particular unit belongs
            idx = find((prodsigns(:,1) .* prodsigns(:,2)) < 0);
            for l = 1 : length(idx)
                poslvl = find(prodsigns(idx(l),:) > 0);
                if poslvl == 1
                    wecod{1,pos}(idx(l),k) = -1 .* (freqtab(tabidx{:}) ./ freqtab(tabidx{poslvl},omit(2)));
                elseif poslvl == 2
                    wecod{1,pos}(idx(l),k) = -1 .* (freqtab(tabidx{:}) ./ freqtab(omit(1),tabidx{poslvl}));
                end
            end
            
        end
        
    end
    
elseif Options.interactions > 2
    
    error('This implementation of weighted effect coding does not handle more than 2-term interactions!');
    
end

% Reconstitutes the full design matrices based on the reduced designs
wecodfull = cell(1,length(wecod));
wecodfull{1} = ones(n,1); % we coding for the grand mean matrix
for i = 2 : length(wecod)
    wecodfull{i} = zeros(n,size(wecod{i},2));
    for j = 1 : size(Ycomb,1)
        idx = ismember(Y, Ycomb(j,:), 'rows');
        wecodfull{i}(idx,:) = repmat(wecod{i}(j,:),[sum(idx),1]);
    end
end

% wecodeempty as an empty cell structure is useful in the main functions
wecodempty = {};
for i = 1 : length(wecodfull)
    wecodempty{i} = zeros(size(wecodfull{i}));
end

end
