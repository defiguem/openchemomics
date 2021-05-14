function [preproc, Options] = xpreproc(data, Options)
%
% The function xpreproc preprocesses data matrices with samples in the
% rows and variables in the columns. The list of the given preprocessing
% techniques available is given below and it is always possible to add new
% ones at will. For each preprocessing technique, the parameters used for
% the calculation are also stored in the output structure. This is
% important, especially when data was divided into calibration and
% validation sets. The parameters calculated with xpreproc for the
% calibration set are then used as inputs in the xpreprop function to
% preprocess the validation set. In that sense, each modification brought
% to the present function must be accompanied with the proper modification
% of the xpreprop function.
%
% Care must be taken in keeping the independence between the calibration
% and validation sets. In fact, some techniques such as the column mean
% centering use the parameters calculated with the calibration set to
% preprocess the validation set. Others, such as SNV work on the individual
% samples and each sample is preprocessed according to its own parameters.
% In the latter case, there is no difference between the calculations made
% with xpreproc or xpreprop.
%
% This function handles preprocessing of one or multiple blocks. Usage must
% be adapted accordingly. For each table, one single or several
% pretreatments can be performed in a row. See usage below for more
% information.
%
% References :
% ============
% Roger, J.-M., Boulet, J.-C., Zeaiter, M., & Rutledge, D. N. (2020).
% 3.01-Pre-processing Methods. In S. Brown, R. Tauler, & B. Walczak (Éds.),
% Comprehensive Chemometrics (Second Edition) (p. 1-75). Elsevier.
% https://doi.org/10.1016/B978-0-12-409547-2.14878-4
%
% Input arguments :
% =================
% data : matrix of samples in the rows and variables in the columns;
%
% Options.prepro.X.type : array of cells containing in each cell the chain
% of pretreatments to apply to each data block;
%
% Options.prepro.X.para : array of cells containing in each cell the
% parameters associated to each pretreatment in Options.prepro.X.type;
%
% Output arguments :
% ==================
% preproc : structure containing the preprocessed data in preproc.data
% along with sub-structures for each of the preprocessing techniques used
% containing the parameters used for calculation.
%
% Usage :
% =======
% Options.prepro.X.type = {'cmeancenter'};
% Options.prepro.X.para = repelem({''},length(Options.prepro.X.type));
% [preprocx] = xpreproc(data, Options);
%
% or
%
% Options.prepro.X.type = {{'autoscale'};{'cmeancenter'}};
% Options.prepro.X.para = repelem({''},length(Options.prepro.X.type));
% [preprocx] = xpreproc(data, Options); % data is a cell array of 2 tables
%
% Related functions :
% ===================
% xpreprop.m (does preprocessing on X testing data)
% ypreproc.m (does preprocessing on Y training data)
% yunprepro.m (does unpreprocessing on Y preprocessed data)
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

% Checks if an array of cells with several tables was provided
if iscell(data) == 1
    ntab = length(data); % number of tables in data input
else
    data = {data}; % place array in cell
    ntab = 1;
end

% Checks if preprocessing methods were chosen for X, otherwise does nothing
if isfield(Options,'prepro') == 0 || isfield(Options.prepro,'X') == 0 || isfield(Options.prepro.X,'type') == 0
    Options.prepro.X.type = repelem({'nothing'},ntab)';
end

% Checks if preprocessing parameters were provided when needed
if isfield(Options.prepro.X,'para') == 0 
    for i = 1 : ntab
        Options.prepro.X.para{1,i} = repelem({''},length(Options.prepro.X.type))';
    end
end

% Checks if preprocessing methods were chosen for X, otherwise does nothing
if isfield(Options.prepro.X,'type') == 1 && ntab ~= length(Options.prepro.X.type)
    error('Preprocessing methods to be applied must be defined for each of the input tables');
end

if ischar(Options.prepro.X.type) == 1 % Cell format needed in the for loop
    Options.prepro.X.type = cellstr(Options.prepro.X.type);
end

%% Main section

% Loops over the block in the input argument data
for t = 1 : ntab
    
    if ischar(Options.prepro.X.type{t}) == 1 % Cell format needed in the for loop
        Options.prepro.X.type{t} = cellstr(Options.prepro.X.type{t});
    end
    
    % Stores the preprocessing steps to be used in the order of application
    preproc(t).preprosteps = Options.prepro.X.type{t};
    
    % Allocates space for storage of preprocessing methods parameters
    preproc(t).prepropara = cell( length(preproc(t).preprosteps), 1 );
    
    % Starts the preprocessing loop with as many iterations as preprocessing
    % techniques specified in Options.prepro.X.type{t}
    for i = 1 : length(Options.prepro.X.type{t})
        
        % Determines the dimensions of the data matrix
        % Needed inside the loop in case some preprocessing techniques modifiy the data size
        [xrow,xcol] = size(data{t});
        
        switch Options.prepro.X.type{t}{i} % Selects one preprocessing technique at a time
            
            case 'nothing'
                % No preprocessing performed
                
            case 'cmeancenter' % Column mean centering (no parameters to optimize)
                cmdata = mean(data{t},1); % Column mean calculation
                data{t} = data{t} - ( ones(xrow,1) * cmdata ); % Substracts from each element in a given column the associated column mean
                
                % Parameter storage
                preproc(t).prepropara{i}.cmeancenter.cmdata = cmdata;
                
            case 'autoscale' % Column autoscaling / variance scaling / column standardization / z-transformation (no parameters to optimize)
                cmdata = mean(data{t},1); % Column mean calculation
                cstddata = std(data{t},0,1); % Column standard deviation calculation
                data{t} = (data{t} - ( ones(xrow,1) * cmdata )) ./ ( ones(xrow,1) * cstddata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation
                
                % Parameters storage
                preproc(t).prepropara{i}.autoscale.cmdata = cmdata;
                preproc(t).prepropara{i}.autoscale.cstddata = cstddata;
                
            case 'tabnorm' % normalization of table to unit variance
                tabnorm = sqrt( sum( sum( data{t}.^2 ) ) ); % frobenius norm
                data{t} = data{t} ./ tabnorm;
                
                % Parameters storage
                preproc(t).prepropara{i}.tabnorm = tabnorm;
                
            otherwise
                error('Preprocessing technique ''%s'' does not exist or is not supported', Options.prepro.X.type{t}{i});
                
        end
        
    end
    
    % The data is preprocessed accross the loop iterations and the final form
    % is stored in the output structure only in the end. However, the
    % parameters calculated throughout the loop iterations are directly stored.
    preproc(t).data = data{t};
    
end

end