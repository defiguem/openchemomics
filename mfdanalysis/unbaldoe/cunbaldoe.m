function [X, Y, sid] = cunbaldoe(X, Y, ns, combid, Options)
%
% The cunbaldoe function unbalances a design of experiments contained in Y.
% The function works even if the input design is already unbalanced but
% it throws a warning if that happens.
%
% Data unbalancing is done in a controlled way. The number of observations 
% to be removed in a given combination of levels in factors must be defined. 
% It means that only one combination of levels is unbalanced, not the others.
%
% Input arguments :
% =================
% X : predictors matrix with observations in the rows and variables in the
% columns;
% 
% Y : design matrix with observations in the rows and factors in the columns
%
% ns : number of observations to remove in a given combination of levels
%
% combid : the identifier of the combination of levels to unbalance
% 
% Options : a structure containing optional user defined parameters
%       Options.rng : seed for random numbers generation
%
% Output arguments :
% ==================
% X : unbalanced input data matrix X
%
% Y : unbalanced design matrix Y
%
% sid : indices of observations removed
%
% Usage :
% =======
% Options.rng = 123456789;
% [Xunbal, Yunbal, sid] = cunbaldoe(Xbal, Ybal, ns, combid, Options)
%
% Related functions :
% ===================
% baldoe.m (balances design matrices)
% runbaldoe.m (unbalances design matrices randomly)
% 
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications :
% ===============
% 
% =========================================================================

%% Fail-safe

[n,p] = size(X); % number of observations and variables
f = size(Y,2); % number of factors

% Checks if the design matrix Y is already unbalanced
for i = 1 : f
    [levels] = histcounts(Y(:,i));
    if length(unique(levels)) ~= 1
        warning(['Factor in column ',num2str(i),'of the design matrix Y is already unbalanced!']);
    end
end

% Checks if the random generator number seed was defined
if isfield(Options,'rng') == 0
    rng('shuffle')
else
    rng(Options.rng); 
end

%% Main section

% Defines all possible unique combinations of factors/levels along with to
% which combination each row of X belongs to
[uni,~,ic] = unique(Y,'rows','stable');

% Calculates how many members per unique combination there are
nuni = zeros(size(uni,1),1); % number of observations per unique combination
for i = 1 : size(uni,1) % number of unique combinations
    nuni(i) = sum(ic == i);
end

rows = find(ismember(Y,uni(combid,:),'rows') == 1);

sid = rows(randperm(length(rows),ns)); 

X(sid,:) = []; % removes observations at sid from X
Y(sid,:) = []; % removes observations at sid from design of experiments Y

end