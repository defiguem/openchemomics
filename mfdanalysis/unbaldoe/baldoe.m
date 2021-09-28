function [Xbal, Ybal, sid] = baldoe(X, Y, Options)
%
% The baldoe function balances a design of experiments contained in Y.
% The function works even if the input design is already balanced but
% it throws a warning if that happens.
%
% Balancing of the design matrix is performed by random undersampling at
% the moment but possible extensions could be considered.
%
% It is possible to define in the options a seed for random number 
% generation. It is useful to make sure that the operation is reproducible.
%
% Input arguments :
% =================
% X : predictors matrix with observations in the rows and variables in the
% columns;
% 
% Y : design matrix with observations in the rows and factors in the columns
% 
% Options : a structure containing optional user defined parameters
%       Options.rng : seed for random numbers generation
%       Options.rsmethod : resampling method to use (only 'rus' available
%       at the moment)
%
% Output arguments :
% ==================
% Xbal : balanced input data matrix X
%
% Ybal : balanced design matrix Y
%
% sid : indices of selected observations to create balanced design
%
% Usage :
% =======
% Options.rng = 123456789;
% [Xbal, Ybal, sid] = baldoe(X, Y, Options)
%
% Related functions :
% ===================
% radecomp.m (resampling ANOVA decomposition)
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

[n,p] = size(X); % number of observations and variables
f = size(Y,2); % number of factors

% Checks if the design matrix Y is already balanced
lgthlev = 0;
for i = 1 : f
    [levels] = histcounts(Y(:,i));
    lgthlev = lgthlev + length(unique(levels));
end

if lgthlev == f
    warning('The design matrix Y is already balanced!');
end

% Checks if the random generator number seed was defined
if isfield(Options,'rng') == 0
    rng('shuffle')
else
    rng(Options.rng); 
end

if isfield(Options,'rsmethod') == 0
    Options.rsmethod = 'rus';
end

%% Main section

% Defines all possible unique combinations of factors/levels along with to
% which combination each row of X belongs to
[uni,~,ic] = unique(Y,'rows','stable');

% Calculates how many observations per unique combination there are
nuni = zeros(size(uni,1),1); % number of observations per unique combination
for i = 1 : size(uni,1) % number of unique combinations
    nuni(i) = sum(ic == i);
end

% Switches case depending on the resampling strategy
switch Options.rsmethod
    
    case 'rus' % random under-sampling
        
        % Number of observations per combination to extract in order to create a balanced dataset
        minuni = min(nuni);
        
        % Allocates space for the balanced matrices and indices
        Xbal = zeros(minuni * size(uni,1),size(X,2));
        Ybal = zeros(minuni * size(uni,1),size(uni,2));
        sid = zeros(minuni * size(uni,1),1);
        
        % Performs a random undersampling of all unique combinations to satisfy miuni
        for i = 1 : size(uni,1)
            
            idx = find(ic == i); % finds observations of the unique combination i
            rp = randperm(length(idx),minuni); % random permutation of the number of member of the unique combination i
            
            % Creates balanced design
            Xbal((i-1) * minuni + 1 : i * minuni, :) = X(idx(rp),:);
            Ybal((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
            sid((i-1) * minuni + 1 : i * minuni) = idx(rp);
            
        end
        
end

end