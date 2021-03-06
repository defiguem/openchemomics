function [radecomp_res] = radecomp(X, Y, Options)
%
% The radecomp function performs a Resampled ANOVA decomposition on the X
% matrix based on the design of experiments in Y. This function is meant to
% be used when the design matrix Y is unbalanced.
%
% Note that several methods for the calculation of the resampled matrix are
% available but most of them lead to the same result. These are just
% different strategies but the most efficient way is to calculate a
% weighted resampled matrix based on the initial dataset.
%
% Input arguments :
% =================
% X : predictors matrix with samples in the rows and variables in the
% columns;
%
% Y : design matrix with samples in the rows and factors in the columns
%
% Options : a structure containing optional user defined parameters
%       Options.rsmethod : the resampling method to use
%           'wavg' : weighted averages (default and fastest way);
%           'exres' : builds exhaustive resamplings in a reduced fashion;
%           'exresfull' : builds the full exhaustive resamplings;
%           'rdmres' : random under-sampling repeated n times
%       if Options.rsmethod = 'rdmres', define :
%           Options.rng : seed for random numbers generation
%           Options.nres : the number of resamplings to perform
%
% Output arguments :
% ==================
% radecomp_res : structure containing results of the Resampling ANOVA
% decomposition (see the adecomp.m function for more information)
%
% Usage :
% =======
% Options.rsmethod = 'wavg'; % default
% [radecomp_res] = radecomp(X, Y, Options); 
%
% Options.rsmethod = 'exres';
% [radecomp_res] = radecomp(X, Y, Options);
%
% Options.rsmethod = 'exresfull';
% [radecomp_res] = radecomp(X, Y, Options);
%
% Options.rsmethod = 'rdmres';
% Options.rng = 123456789;
% Options.nres = 100;
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

% Checks if the resampling method was chosen
if isfield(Options,'rsmethod') == 0
    Options.rsmethod = 'wavg';
end

% Checks if the number of resamplings to perform was defined
if strcmp(Options.rsmethod,'rdmres') == 1 && isfield(Options,'nres') == 0
    Options.nres = 100;
    nres = Options.nres;
elseif strcmp(Options.rsmethod,'rdmres') == 1 && isfield(Options,'nres') == 1
    nres = Options.nres;
end

% Checks if the random generator number seed was defined
if isfield(Options,'rng') == 0
    Options.rng = 123456789;
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

% Number of samples per combination to extract in order to create a balanced dataset
minuni = min(nuni);

switch Options.rsmethod
    
    case 'wavg' % weighted average (faster)
        
        Xbal = zeros(minuni * length(nuni),size(X,2));
        Ybal = zeros(minuni * size(uni,1),size(Y,2));
        for i = 1 : size(uni,1)
            Xbal((i-1) * minuni + 1 : i * minuni,:) = repmat(mean(X(ic == i,:),1),[minuni,1]);
            Ybal((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
        end
        
    case 'exres' % exhaustive resampling
        
        % Calculates all possible combinations of observations for each
        % possible combination of levels
        nres = cell(1,length(nuni));
        maxnres = 0;
        for i = 1 : length(nuni)
            nres{i} = nchoosek(1:nuni(i),minuni);
            nres{i} = nres{i}(randperm(size(nres{i},1)),:);
            if size(nres{i},1) > maxnres
                maxnres = size(nres{i},1);
            end
        end
        
        % Allocates space for the output results
        Xbal = zeros(minuni * size(uni,1),size(X,2));
        Ybal = zeros(minuni * size(uni,1),size(Y,2));
        
        % Loops over the unique combinations of levels
        for i = 1 : length(nres)
            
            idx = find(ic == i); % finds observations of the unique combination i
            
            % Calculates a resampled matrix for each combination of samples
            % (the resampled matrices are cumulated in a 2D matrix and averaged
            % later for memory space concerns in case size(nres{i},1) is large)
            Xres = zeros(minuni,size(X,2));
            for j = 1 : size(nres{i},1)
                nres{i}(j,:) = idx(nres{i}(j,:));
                Xres = Xres + X(nres{i}(j,:),:);
            end
            
            % Stores the resampled mean matrix
            Xbal((i-1) * minuni + 1 : i * minuni,:) = Xres ./ size(nres{i},1);
            
            % Stores the balanced design of experiments
            Ybal((i-1) * minuni + 1 : i * minuni, :) = repmat(uni(i,:),[minuni,1]);
            
        end
        
        % Stores the number of resamplings performed
        Options.nres = nres;
        
    case 'exresfull' % exhaustive resampling
        
        % Calculates all possible combinations of observations for each
        % possible combination of levels
        nres = cell(1,length(nuni));
        maxnres = 0;
        for i = 1 : length(nuni)
            nres{i} = nchoosek(1:nuni(i),minuni);
            nres{i} = nres{i}(randperm(size(nres{i},1)),:);
            if size(nres{i},1) > maxnres
                maxnres = size(nres{i},1);
            end
            
            for j = 1 : size(nres{i},1)
                idx = find(ic == i); % finds observations of the unique combination i
                nres{i}(j,:) = idx(nres{i}(j,:));
            end
            
        end
        
        C = cellfun(@(x) (1:size(x,1)),nres,'UniformOutput',false);
        C = combvec(C{:});

        % Allocates space for the output results
        Xbal = [];
        Ybal = [];
            
        % Loops over the unique combinations of levels
        for i = 1 : size(C,2)
            
            for j = 1 : size(C,1)
                Xbal = [Xbal;X(nres{j}(C(j,i),:),:)];
                Ybal = [Ybal;Y(nres{j}(C(j,i),:),:)];
            end
            
        end
        
        % Stores the number of resamplings performed
        Options.nres = nres;
    
    case 'rdmres' % n iterations of resampling (slower)
        
        % Allocates space for the output results
        Xbalall = zeros(minuni * size(uni,1),size(X,2),nres);
        sids = zeros(minuni * size(uni,1),nres);
        
        % Loops over the resampling iterations
        for i = 1 : nres
            
            % Balances the design
            Options.rng = Options.rng + 1;
            [Xbal, Ybal, sid] = baldoe(X, Y, Options);
            
            % Stores the resampled matrix
            Xbalall(:,:,i) = Xbal;
            
            % Stores the indices of removed samples
            sids(:,i) = sid;
            
        end
        
        % Calculates a mean resampling matrix
        Xbal = nanmean(Xbalall,3);
        
end

% Performs ANOVA decomposition on the mean resampling matrix
[radecomp_res] = adecomp(Xbal, Ybal, Options);

end