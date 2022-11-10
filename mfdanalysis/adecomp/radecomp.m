function [radecomp_res, Xbal, Ybal] = radecomp(X, Y, Options)
%
% The radecomp function performs a Rebalanced ANOVA decomposition on the X
% matrix based on the design of experiments in Y. This function is meant to
% be used when the design matrix Y is unbalanced. It can also be used if Y
% is balanced but in this case, Xbal and Ybal will be equal to the inputs X
% and Y.
%
% Note that several methods for the calculation of the rebalanced matrix 
% are available but most of them lead to the same result. 
%
% Reference :
% ===========
% de Figueiredo, M., Giannoukos, S., Rudaz, S., Zenobi, R., & Boccard, J. 
% (s. d.). Efficiently handling high-dimensional data from multifactorial 
% designs with unequal group sizes using Rebalanced ASCA (RASCA). 
% Journal of Chemometrics, n/a(n/a), e3401. https://doi.org/10.1002/cem.3401
%
% Input arguments :
% =================
% X : predictors matrix with samples in the rows and variables in the
% columns;
%
% Y : design matrix with samples in the rows and factors in the columns
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
% Usage :
% =======
% Options.rebmethod = 'overavg'; % default
% [radecomp_res] = radecomp(X, Y, Options);
%
% Options.rebmethod = 'underavg'; 
% [radecomp_res] = radecomp(X, Y, Options); 
%
% Related functions :
% ===================
% adecomp.m (performs ANOVA decomposition)
% sumcoding.m (sum coding of the design matrices for GLM)
% wecoding.m (weighted-effect coding of the design matrices for GLM)
% baldoe.m (balances design matrices)
% unbaldoe.m (unbalances design matrices)
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

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if the rebalancing method was chosen
if isfield(Options,'rebmethod') == 0
    Options.rebmethod = 'overavg';
end

%% Main section

% Defines all possible unique combinations of factors/levels along with to
% which combination each row of X belongs to
[uni,~,ic] = unique(Y,'rows','stable');

% Calculates how many member per unique combination there are
nuni = zeros(size(uni,1),1); % number of samples per unique combination
for i = 1 : size(uni,1) % number of unique combinations
    nuni(i) = sum(ic == i);
end

switch Options.rebmethod
    
    case 'underavg' % averages 

        % Number of samples per combination to extract in order to create a balanced dataset
        minuni = min(nuni);
        
        % Allocates space for the balanced design and data
        Xbal = zeros(minuni * length(nuni),size(X,2));
        Ybal = zeros(minuni * size(uni,1),size(Y,2));

        % Loops over each unique combination of levels
        for i = 1 : size(uni,1)
            Xbal((i-1) * minuni + 1 : i * minuni,:) = repmat(mean(X(ic == i,:),1),[minuni,1]);
            Ybal((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
        end
           
    case 'overavg'
        
        % Number of samples per combination to extract in order to create a balanced dataset
        minuni = max(nuni);

        % Allocates space for the balanced design and data
        Xbal = zeros(minuni * length(nuni),size(X,2));
        Ybal = zeros(minuni * size(uni,1),size(Y,2));

        % Loops over each unique combination of levels
        for i = 1 : size(uni,1)
            idx = (ic == i); % finds all replication of combinatio  i
            if sum(idx) == minuni % if already enough replicates
                Xbal((i-1) * minuni + 1 : i * minuni,:) = X(idx,:);
            else % if missing replicates
                Xbal((i-1) * minuni + 1 : i * minuni,:) = [X(idx,:);repmat(mean(X(idx,:),1),[minuni-sum(idx),1])];
            end
            Ybal((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
        end
           
end

% Performs ANOVA decomposition on the rebalanced matrix
[radecomp_res] = adecomp(Xbal, Ybal, Options);

end