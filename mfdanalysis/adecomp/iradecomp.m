function [radecomp_res, Xs, Ys] = iradecomp(X, Yrep, Y, Options)
% 
% The iradecomp function performs a Rebalanced ANOVA decomposition on the 
% tables in X based on the designs of experiments in Y. This function is 
% meant to be used when the design matrix Y is unbalanced. It can also be 
% used if Y is balanced but in this case, the function adecomp can be used.
%
% This function handles the rebalancing in two-stages:(1) the tables are
% homogenized as to make sure that the same row number in all tables
% relate to the same replicate. (2) the design is balanced if needed.
%
% Reference :
% ===========
%
% Input arguments :
% =================
% X : array of cells with tables in which predictors matrix with samples 
% in the rows and variables in the columns;
%
% Yrep : array of cells with the identity of each replicate as an integer;
%
% Y : array of cells with the design matrix of each table with samples in 
% the rows and factors in the columns;
%
% Options : a structure containing optional user defined parameters
%       Options.rebmethod : the rebalancing method to use
%           'underavg' : calculates averages for each unique combination of
%                        levels from the available replicates.
%           'overavg' : finds unique combinations of levels with missing
%                       replicates and replaces them with the average of 
%                       all measured replicates associates with that 
%                       combination of levels. This method has the
%                       advantage of retaining data's natural veriation and
%                       should be prefered (default).
%
% Output arguments :
% ==================
% radecomp_res : structure containing results of the Rebalanced ANOVA
% decomposition (see the adecomp.m function for more information)
%
% Xs : the same as input X aligned wrt potentially missing replicates;
%
% Ys : the same as input Y aligned wrt potentially missing replicates;
%
% Usage :
% =======
% Options.rsmethod = 'overavg'; % default
% [radecomp_res] = iradecomp(X, Y, Options);
%
% Options.rsmethod = 'underavg'; 
% [radecomp_res] = iradecomp(X, Y, Options); 
%
% Related functions :
% ===================
% adecomp.m (performs ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications :
% ===============
%
% =========================================================================

%% Fail-safe section

ntab = length(X); % number of tables

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if the rebalancing method was chosen
if isfield(Options,'rebmethod') == 0
    Options.rebmethod = 'overavg';
end

% Checks that X, Y and Yrep have the same number of cells
if length(X) ~= length(Yrep) || length(X) ~= length(Y)
    error('The number of cells in X, Y and Yrep must be equal. Please check.');
end

% Checks that the number of rows in each cell of X, Y and Yrep are equal,
% otherwise it throws an error.
for i = 1 : ntab
    if size(X{i},1) ~= length(Yrep{i}) || size(X{i},1) ~= size(Y{i},1)
        error('The number of rows in each cell of X, Y and Yrep must be equal. Please check.');
    end
end

% Nota that we allow the tables in X to have different numbers of 
% replicates, but we impose the design Y and replicates Yrep to have the
% same number of rows as the associated cell in X.

%% Equilibrates the designs as to have each replicate appear in each of the tables

% Usually for multiblock methods, the same row in each table relates to the
% same replicate. If for some reason this is not the case, this section
% aims at aligning the rows in the tables wrt to the information in Yrep.
% This is the first step towards balancing the design.

% Identifies the highest number of replicates
maxrep = max(cellfun(@(x) max(x), Yrep));
refrep = (1 : maxrep)'; % reference vector for the replicates

% Loops over each table
Xs = cell(1, ntab);
Ys = cell(1, ntab);
nans = zeros(length(refrep),ntab);
for i = 1 : ntab

    % Allocates space with NaNs
    Xs{i} = nan(length(refrep),size(X{i},2));
    Ys{i} = nan(length(refrep),size(Y{i},2));

    % Fills the NaN arrays with the measured values
    idx = ismember(refrep, Yrep{i});
    Xs{i}(idx,:) = X{i};
    Ys{i}(idx,:) = Y{i};

    % Stores the position of the missing replicates in each table
    nans(~idx,i) = 1; % to be placed in output

end

% Creates one global Y to defines the design
Yall = mean(cat(3,Ys{:}),3,'omitnan');

% Note that we equilibrated the design between the different tables, but it
% doesn't mean that the exprimental design is balanced. There is a shared
% design between all tables defined now by Yall

%% Fills values for missing replicates missing 

% Finds the position of the replicates that were missing
[r,c] = find(nans == 1);

% Loops over missing replicates and assigns mean of associated unique combination of levels
for i = 1 : length(r)
    idx = ismember(Yall, Yall(r(i),:),'rows'); % finds all replicates with the same unique combination of levels
    Xs{c(i)}(r(i),:) = mean(Xs{c(i)}(idx,:),'omitnan'); % assigns the mean of these replicates, ignoring NaNs
end

for i = 1 : length(Y)
    Ys{i} = Yall;
end

%%  Finds the rows missing in all tables

idx = find(sum(nans,2) == size(nans,2));
if isempty(idx) == 0
    for i = 1 : ntab
        Xs{i}(idx,:) = [];
        Ys{i}(idx,:) = [];
    end

    Yall(idx,:) = [];
end

%% Balances the design

% Defines all possible unique combinations of factors/levels along with to
% which combination each row of X belongs to
[uni,~,ic] = unique(Yall,'rows','stable');

% Calculates how many member per unique combination there 
nuni = zeros(size(uni,1),1); % number of samples per unique combination
for i = 1 : size(uni,1) % number of unique combinations
    nuni(i) = sum(ic == i);
end

% Allocates space for the rebalanced design of all tables 
Xbal = cell(1, ntab);
Ybal = cell(1, ntab);

switch Options.rebmethod
    
    case 'underavg' % averages 

        % Number of samples per combination to extract in order to create a balanced dataset
        minuni = min(nuni);
        
        % Loops over each table
        for j = 1 : ntab

            % Allocates space for the balanced design and data
            Xbal{j} = zeros(minuni * length(nuni),size(Xs{j},2));
            Ybal{j} = zeros(minuni * size(uni,1),size(Ys{j},2));

            % Loops over each unique combination of levels
            for i = 1 : size(uni,1)
                Xbal{j}((i-1) * minuni + 1 : i * minuni,:) = repmat(mean(Xs{j}(ic == i,:),1),[minuni,1]);
                Ybal{j}((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
            end

        end
           
    case 'overavg'
        
        % Number of samples per combination to extract in order to create a balanced dataset
        minuni = max(nuni);

        % Loops over each table
        for j = 1 : ntab

            % Allocates space for the balanced design and data
            Xbal{j} = zeros(minuni * length(nuni),size(Xs{j},2));
            Ybal{j} = zeros(minuni * size(uni,1),size(Ys{j},2));

            % Loops over each unique combination of levels
            for i = 1 : size(uni,1)
                idx = (ic == i); % finds all replication of combinatio  i
                if sum(idx) == minuni % if already enough replicates
                    Xbal{j}((i-1) * minuni + 1 : i * minuni,:) = Xs{j}(idx,:);
                else % if missing replicates
                    Xbal{j}((i-1) * minuni + 1 : i * minuni,:) = [Xs{j}(idx,:);repmat(mean(Xs{j}(idx,:),1),[minuni-sum(idx),1])];
                end
                Ybal{j}((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
            end
        end
           
end

% Performs ANOVA decomposition on the rebalanced matrix of each table
for i = 1 : ntab
    [radecomp_res(i)] = adecomp(Xbal{i}, Ybal{i}, Options);
end

% Adds the information on the NaNs to the structure array
nans = mat2cell(nans,size(nans,1),ones(1,size(nans,2)));
[radecomp_res.nans] = nans{:};


end