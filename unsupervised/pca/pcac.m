function [pca_model] = pcac(X, PCs, Options)
%
% The pcac function performs Principal Components Analysis (PCA) using
% the Singular Value Decomposition (SVD) by default and possibly other 
% methods depending on the choices made by the user.
% 
% Note that no preprocessing is applied by the function and any required
% preprocessing should be performed before.
% 
% Input argument:
% ===============
% X : data matrix of samples placed in the rows and variables placed in the columns;
% 
% PCs : number of principal components to extract
% 
% Output argument:
% ================
% pca_model: structure containing results of the PCA 
%   pca_model.scores : the scores of the PCA
%   pca_model.loadings : the loadings (eigenvectors) of the PCA
%   pca_model.evals : the eigenvalues of the eigenvectors
%   pca_model.varexp : percentage of explained variances the Principal
%   Components (PCs or eigenvectors)
%   pca_model.cumvarexp : Cumulative percentage of explained variance by the
%   PCs
% 
% USAGE :
% =======
% Options.pcatype = 'svd'; % or 'nipals' or 'evd'
% [pca_model] = pcac(X, PCs, Options);
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

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if pca method was defined
if isfield(Options,'pcatype')==0
    Options.pcatype = 'svd';
end

[n,p] = size(X); % size of X

if PCs > min([n,p]) % maximum number of PCs to extract
    PCs = min([n,p]); % number of PCs to extract
end

%% Main section

switch Options.pcatype % cases for performing PCA
    
    case 'svd'
        
        % SVD function with the output parameters U,S and V
        % 'econ' avoids unnecessary zeros provinding compact form matrices
        [U,S,V] = svd(X,'econ');
        
        % Calculates the objects' scores in the new referential system
        T = U * S;
        
        % Calculates the eigenvalues by squarring the diagonal elements of
        % the S matrix containing the singular values
        evals = S.^2;
        
    case 'nipals' % Nonlinear Iterative Partial Least Squares
        
        PCs = min([n,p]);
        
        itlim = 1000; % arbitrary nb of iterations for each PC extraction
        tol = 1e-6; % arbitrary tolerance nb between two while loops
        
        E = X; % matirx of residuals set to be X at first

        evals = zeros(PCs,PCs); % allocates space for the eigenvalues
        T = zeros(n,PCs); % allocates space for the scores
        V = zeros(p,PCs); % allocates space for the loadings
        
        for i = 1 : PCs
            
            it = 1;
            told = E(:,1); % initialization of t by first column in X
            
            while 1 % while loop will only stop if break
                
                t = told; % useful for convergence verification
                p = (E'*t)/(t'*t); % Projection of E onto t to find the corresponding loadings p
                p = p/norm(p); % normalizes loadings p to have length 1
                t = (E*p)/(p'*p); % projection of E onto p to find the corresponding scores t
                
                if (norm(t-told) < tol) || (it == itlim)
                    break;
                end
                
                told=t;
                it=it+1;
                
            end
            
            T(:,i) = t; % stores the scores iteratively calculated
            V(:,i) = p; % stores the loadings iteratively calculates
            evals(i,i) = trace(((t*p')'*(t*p'))); % calculates the eigenvalues as the trace of the covariance matrix
            E = E - t*p'; % calculates the residual matrix
            
        end
        
    case 'evd'
        
        % Performs Eigenvalue Decomposition
        [V,evals] = eig(X' * X);
        
        % Sorts the eigenvalue indecreasing order
        [evals,idx] = sort(diag(evals),'descend');
        evals = diag(evals);
        
        % Sorts the eigenvectors (loadings)
        V = V(:,idx);
        
        % Calculates the scores
        T = X * V;
        
end

% Creates the ouput structure containing the output parameters
pca_model.PCs = PCs; % number of principal components
pca_model.scores = T(:,1:PCs); % scores
pca_model.loadings = V(:,1:PCs); % loadings
pca_model.evals = diag( evals(1:PCs,1:PCs) ) ./ (n-1); % Eigenvalues
pca_model.varexp = 100 * ( pca_model.evals ./ sum(diag(evals) ./ (n-1)) ); % Calculates the explained variance (varexp) by the PCs
pca_model.cumvarexp = cumsum(pca_model.varexp); % Calculates the cumulative sum of the explained variance by the PCs

pca_model.Qcont = X - ( T(:,1:PCs) * V(:,1:PCs)' ); % residuals for the X matrix (contributions)
pca_model.Qresobj = sum(pca_model.Qcont.^2,2); % Objects' redisuals

S = diag(evals) ./ (n-1); 

% Contributions of Hotteling's T2 / Leverage
T = pca_model.scores;
pca_model.Tcont = T*diag(sqrt(1./S(1:PCs)))*pca_model.loadings';
for i=1:n
    pca_model.Thotobj(i,1) = T(i,:)*pinv(diag(S(1:PCs)))*T(i,:)';
end
pca_model.Tlimobj = ((PCs*(n-1))/((n-PCs)))*finv(0.95,PCs,n-PCs); % The UCL (Upper Control Limit) is here fixed at 0.95 or 95%

% Calculation of residuals' confidence limit according to Jackson-Mudholkar
if PCs >= 2
    S = S(S>0); theta1=sum(S(PCs+1:end),1); theta2=sum(S(PCs+1:end).^2,1); theta3=sum(S(PCs+1:end).^3,1); % Calculation of the theta parameters
    h0 = 1-(2*theta1*theta3)/(3*theta2^2); c_alpha = norminv(0.95,0,1);
    pca_model.Qlimobj=theta1*((c_alpha*sqrt(2*theta2*h0^2)/theta1)+(theta2*h0*(h0-1)/theta1^2)+1)^(1/h0);
else
    % Calculation of confidence limit for X residuals
    DF = 2 * ( mean(pca_model.Qresobj) / std(pca_model.Qresobj) )^2; % Simpler alternative for Qlim : pca_model.Qlimobj=prctile(pca_model.Qresobj,0.95*100);
    pca_model.Qlimobj = chi2inv(0.95,DF) * mean(pca_model.Qresobj) / DF; % According to http://mdatools.com/mdatools/residuals-and-critical-limits.html
end

pca_model.Options=Options;

end