function [X, Y, sid] = runbaldoe(X, Y, Options)
%
% The runbaldoe function unbalances a design of experiments contained in Y.
% The function works even if the input design is already unbalanced but
% it throws a warning if that happens. Data unbalancing is done randomly
% according to a given percentage of the data to remove.
%
% Input arguments :
% =================
% X : predictors matrix with samples in the rows and variables in the
% columns;
% 
% Y : design matrix with samples in the rows and factors in the columns
% 
% Options : a structure containing optional user defined parameters
%       Options.percent : percentage of samples to remove randomly;
%       Options.rng : seed for random numbers generation
%
% Output arguments :
% ==================
% X : unbalanced input data matrix X
%
% Y : unbalanced design matrix Y
%
% sid : indices of samples removed
%
% Usage :
% =======
% Options.rng = 123456789;
% Options.percent = 2;
% [Xunbal, Yunbal, sid] = runbaldoe(X, Y, Options)
%
% Related functions :
% ===================
% baldoe.m (balances design matrices)
% cunbaldoe.m (unbalances design matrices in a controlled way)
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

[n,p] = size(X); % number of samples and variables
f = size(Y,2); % number of factors

% Checks if the design matrix Y is already unbalanced
for i = 1 : f
    [levels] = histcounts(Y(:,i));
    if length(unique(levels)) ~= 1
        warning(['Factor in column ',num2str(i),'of the design matrix Y is already unbalanced!']);
    end
end

% Checks if the percentage of samples to remove was defined
if isfield(Options,'percent') == 0
    Options.percent = 10; 
end

% Checks if the random generator number seed was defined
if isfield(Options,'rng') == 0
    rng('shuffle')
else
    rng(Options.rng); 
end

%% Main section

sid = randperm(n,floor(Options.percent / 100 * n)); % random indices of samples to remove

X(sid,:) = []; % removes samples at sid from X
Y(sid,:) = []; % removes samples at sid from design of experiments Y

end