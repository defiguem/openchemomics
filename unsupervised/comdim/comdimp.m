function [comdim_res] = comdimp(Xs, comdim_model)
%
% The comdimp function takes a calibrated ComDim model to project new
% samples into it. The model input is the output of the comdimc.m function.
%
% Note that this function does not perform preprocessing of the tables 
% according to the calibration preprocessing. The user is invited to see 
% function comdimcp.m to see how preprocessing is handled. All tables 
% should be normed to unit variance to avoid that individual tables with 
% higher variance have too much influence during the extraction of CCs.
%
% Input arguments :
% =================
% Xs : array of cells with as many cells as there are X tables. The number
%   of tables must be consistent with the one used to build the model
%
% comdim_model : ComDim calibrated model with comdimc.m function 
%
% Output arguments :
% ==================
% comdim_res : structure of fields containing projections into ComDim model
%   comdim_res.Q : projected global scores
%   comdim_res.Qloc : projected local scores
% 
% Usage :
% =======
% [comdim_res] = comdimp(Xs, comdim_model);
%
% Related function :
% ==================
% comdimc.m (calibrates ComDim model)
% comdimcp.m (calibrates and predicts samples, preprocessing included)
% 
% Author :
% Miguel de Figueiredo
% @ : miguel.defigueiredo@org.chem.ethz.ch
% 
% Modifications:
% ==============
%
% =========================================================================

%% Fail-safe section

ntabc = size(comdim_model.saliences, 1); % number of calibration tables
ntabp = length(Xs); % number of prediction tables

% Checks if the number of tables for prediction is the same than calibration
if ntabc ~= ntabp 
    error('Error : the number of tables must be the same for calibration and prediction');
end

% Checks the size of the tables in Xs, which need to be the same
for i = 1 : ntabp; xrow(i) = size(Xs{1,i}, 1); end
if numel(unique(xrow)) ~= 1
    error('Error in comdimp.m : data tables in Xs must all have the same number of rows');
else
    xrow = unique(xrow); % defines the number of samples per table
end

%% Main section

CCs = comdim_model.CCs; % extracts the number of CCs from comdim input model
P = comdim_model.P; % extracts global loadings from comdim input model
Ploc = comdim_model.Ploc; % extracts local loadings from comdim input model

Q = zeros(xrow,CCs); % allocates space for prediction global scores
Qloc = cell(1,ntabc); % allocates space for predictionn local scores
Xw = cell(1,ntabc); % allocates space for the saliences-weighted input Xs

for j = 1 : CCs % main loop to extract CCs
    
    % Loops over tables and weights Xs with saliences
    for i = 1 : ntabc 
        Xw{1,i} = comdim_model.saliences(i,j) * Xs{1,i}; 
    end
    
    % Calculates and standardizes the global scores
    Q(:,j) = cat(2,Xw{:}) * P(:,j) * pinv(P(:,j)' * P(:,j));
    Q(:,j) = Q(:,j) / sqrt(Q(:,j)' * Q(:,j)); 
    
    % Calculates local scores with the unweighted Xs
    % Then, deflates the Xs after each CC extraction
    for i = 1 : ntabc
        Qloc{1,i}(:,j) = Xs{1,i} * Ploc{1,i}(:,j) * pinv(Ploc{1,i}(:,j)' * Ploc{1,i}(:,j)); 
        Xs{1,i} = Xs{1,i} - ( Q(:,j) * Ploc{1,i}(:,j)' ); 
    end

end

%% Output structure

comdim_res.Q=Q; % global scores
comdim_res.Qloc=Qloc; % local scores

end