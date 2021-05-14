function [comdim_model] = comdimc(Xs, CCs, Options)
%
% The comdimc function calibrates a ComDim model for the tables stored in
% the Xs array of cells. Each table in Xs must have the same number of rows
% but can have different numbers of variables.
%
% ComDim finds dimensions of greater dispersion common to all tables. These
% dimensions are called Common Components (CCs).
%
% Note that this function does not perform preprocessing of the tables. The
% user is invited to see function comdimcp.m to see how preprocessing is 
% handled. All tables should be normed to unit variance to avoid that
% individual tables with higher variance have too much influence during the
% extraction of CCs.
%
% References :
% ============
% Qannari, E. M., Wakeling, I., Courcoux, P., & MacFie, H. J. H. (2000). 
% Defining the underlying sensory dimensions. Food Quality and Preference, 
% 11(1-2), 151-154. https://doi.org/10.1016/S0950-3293(99)00069-5
%
% Mazerolles, G., Hanafi, M., Dufour, E., Bertrand, D., & Qannari, E. M. 
% (2006). Common components and specific weights analysis : A chemometric 
% method for dealing with complexity of food products. Chemometrics and 
% Intelligent Laboratory Systems, 81(1), 41-49. 
% https://doi.org/10.1016/j.chemolab.2005.09.004
%
% Cariou, V., Jouan-Rimbaud Bouveresse, D., Qannari, E. M., & 
% Rutledge, D. N. (2019). ComDim Methods for the Analysis of Multiblock 
% Data in a Data Fusion Perspective. In Data Handling in Science and 
% Technology (Vol. 31, pp. 179-204). 
% https://doi.org/10.1016/B978-0-444-63984-4.00007-7
%
% Acknowledgments :
% =================
% Acknowledgments to Prof. Emeritus Douglas Neil Rutledge for the
% introduction to ComDim methods and their implementations.
%
% Input arguments :
% =================
% Xs : cell arrays with as many cells as there are tables in Xs
%
% CCs : number of common components to extract
%
% Options : options to be defined by the user, otherwise sets to default
%       Options.convcrit : convergence criterion to extract each CC
%
% Output arguments :
% ==================
% comdim_model : structure of fields containing the ComDim model
%   comdim_model.CCs : number of CCs extracted
%   comdim_model.Q : global scores
%   comdim_model.P : global loadings
%   comdim_model.Qloc : local scores
%   comdim_model.Ploc : local loadings
%   comdim_model.saliences : weight of the tables in the CCs construction
%   comdim_model.varexp : variance explained by CCs
%   comdim_model.Options : options used
% 
% Usage :
% =======
% [comdim_model] = comdimc(Xs, CCs, Options);
%
% Related function :
% ==================
% comdimp.m (projects samples into the space defined by the ComDim model)
% comdimcp.m (calibrates and predicts samples, preprocessing included)
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

ntab = size(Xs,2); % defines the number of tables
n = size(Xs{1},1);

% Checks the size of the tables in Xs, which need to be the same
for i = 1 : ntab; xrow(i) = size(Xs{1,i}, 1); end
if numel(unique(xrow)) ~= 1
    error('Error in comdimc.m : data tables in Xs must all have the same number of rows');
else
    xrow = unique(xrow); % defines the number of samples per table
end

% Checks if a Options structure exists, if not creates it
if exist('Options','var') == 0
    Options = struct;
end

% Checks if a convergence criterion was defined, otherwise sets to default
if isfield(Options,'convcrit') == 0
    convcrit = 1e-10; % convergence criterion to extract CCs
    Options.convcrit = 1e-10; % saves parameter to options
end
    
%% Main section

Xsvo = Xs; % creates a copy of the input tables because Xs will be deflated
saliences = zeros(ntab, CCs); % allocates space for saliences
Q = zeros(xrow, CCs); % allocates space for global scores
% Allocates space for other ComDim parameters
Qloc = cell(1,ntab); % local scores
Ploc = cell(1,ntab); % local loadings

% Calculates the total sum of squares for all blocks Xs 
ssx = zeros(ntab, 1);
for i = 1 : ntab
    ssx(i) = sum(Xs{i}(:).^2);
end

% Calculates associate matrices
Wb = cell(1,ntab);
for i = 1 : ntab
    Wb{i} = Xs{1,i} * Xs{1,i}'; % calculates (X*X')
end

% Allocates space for the explained variance
varexp = struct;
varexp.Xs = zeros(CCs, ntab);

% Loops over the common components to extract
for i = 1 : CCs
    
    % Saliences of tables all equal to 1 at first
    lambda = ones(1,ntab); 
    
    % Initialization of the global scores, normalized to unit length
    q = -1 + 2.*rand(n,1); q = q/sqrt(q'*q);
    
    % Initializes old global scores to enter while loop
    qold=100;
    
    % Runs until stability of global scores between successive iterations
    while norm(qold-q) > convcrit
        
        qold = q; % updates old global scores
                
        % Loops over tables
        Ws = zeros(n,n);
        for j = 1 : ntab
            Ws = Ws + lambda(j) * Wb{j};
        end
        
        [u,sv,v] = svd(Ws,'econ');
        q = u(:,1);
        
        % Correction for sign if needed
        if abs( min(q) ) > abs( max(q) )
            q = -q;
        end

        for j = 1 : ntab
            lambda(j) = q' * Wb{j} * q; % calculates lambda for each table
        end

    end

    saliences(:,i) = lambda; % stores the salience for the CC_i
    Q(:,i) = q; % stores the global scores for the CC_i
    
    for j = 1 : ntab
        Ploc{j}(:,i) = q' * Xs{j};
        Qloc{j}(:,i) = Xs{1,j} * Ploc{j}(:,i) * pinv(Ploc{j}(:,i)' * Ploc{j}(:,i)); 
        Xhat = q * Ploc{j}(:,i)';
        varexp.Xs(i,j) = sum( Xhat(:).^2 );
        Xs{j} = Xs{j} - Xhat;
        Wb{j} = Xs{1,j} * Xs{1,j}';
    end

    % Calculates global explained variance by common component for all Xs
    varexp.Xsglobal(i) = sum(varexp.Xs(i,:)) ./ sum(ssx) .* 100;
    for j = 1 : ntab
        varexp.Xs(i,j) = varexp.Xs(i,j) ./ ssx(j) .* 100;
    end

end

% Calculation of the scaled global loadings
P = cat(2, Xsvo{:})' * Q; % scores are weighted by saliences here

%% Output structure 

comdim_model.CCs = CCs; % number of common components extracted
comdim_model.Q = Q; % global scores
comdim_model.P = P; % global loadings
comdim_model.Qloc = Qloc; % local scores
comdim_model.Ploc = Ploc; % local loadings
comdim_model.saliences = saliences; % saliences of tables
comdim_model.varexp = varexp; % explained variance
comdim_model.Options = Options; % options used

end
