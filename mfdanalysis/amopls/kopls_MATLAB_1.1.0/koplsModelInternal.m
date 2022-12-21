function outval=koplsModelInternal(Xtr,Ytr,A,oax,nrcv,cvType,preProcK,preProcY,cvFrac,modelType,kernelType,Xnew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to make a K-OPLS model for nested cross validation purposes, 
% with returning output (1-Q2) for regression and  (1-sensitivity) for 
% discriminant analysis target function, or (1-area under ROC curve) for 
% 'daAUC' discriminant analysis. These outputs are used as minimisation 
% criteria in grid search and simulated annealing optimisation of the 
% kernel parameter setting, see also 'koplsCVopt'.
%
% Use:
% [outval]=koplsModelInternal(Xtr,Ytr,A,oax,nrcv,cvType,preProcK,preProcY,
%                             cvFrac,modelType,kernelType,Xnew)
%
%** INPUT
% Xtr = The kernel X training set.
% Ytr = The response matrix Y training set.
% A = The number of Y-predictive components (integer).
% oax = The number of Y-orthogonal components (integer).
% nrcv = Number of cross-validation rounds (integer).
% cvType = Type of cross-validation. Either 'nfold' for n-fold
%   cross-validation, 'mccv' for Monte Carlo CV or 'mccvb' for
%   Monte Carlo class-balanced CV. 
% preProcK = Pre-processing settings for the kernel matrix.
%   Either 'mc' for mean-centering or 'no' for no pre-processing. 
% preProcY = Pre-processing parameter for Y. Either 'mc' for
%   mean-centering, 'uv' for mc + scaling to unit-variance,
%   'pa' for mc + Pareto-scaling or 'no' for no scaling.
% cvFrac = Fraction of observations in the training set during
%   cross-validation. Only applicable for 'mccv' or 'mccvb'
%   cross-validation (see 'cvType').
% modelType = 're' for regression, 'da' and 'daAUC' for discriminant 
%   analysis, if 'da', the mean sensitivity is optimised, for 'daAUC' the 
%   area under the receiver operating characteristic curve is optimised 
%   (only for two-class problems). 
% kernelType = kernel type, e.g. 'g' for Gaussian or 'p' for polynomial
%   When using 'p' make sure the kernel parameter is not too big
%   to prevent the polynomial from exploding, resulting in SVD errors.
% Xnew = the kernel parameter used for evaluation.
%
% ** OUTPUT
% outval
%   for discriminant analysis ('da')    = 1- sensitiviy
%   for discriminant analysis ('daAUC') = 1-area under ROC curve.
%   for regression analysis ('re)       = 1-Q^2
%    These values are calculated from a kernel OPLS model using
%    cross-validation without kernel optimisation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Max Bylesjö, Umeå University and Judy Fonville and Mattias 
% Rantalainen, Imperial College.
%   
% Copyright (c) 2007-2010 Max Bylesjö, Judy Fonville and Mattias
% Rantalainen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file is part of the K-OPLS package.
%
% The K-OPLS package is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% The K-OPLS package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelCV=koplsCVopt(Xtr,Ytr,A,oax,modelType,{'nrcvouter',nrcv,'cvType',cvType,'preProcK',preProcK,'preProcY',preProcY,'cvFrac',...
    cvFrac,'verbose',0,'kernelType',kernelType,'kernelParams',Xnew,'opt','no'});
if strcmp(modelType,'re')
    outval=1-modelCV.cv.Q2Yhat(end); 
elseif strcmp(modelType,'da')
    outval=1-modelCV.da.meanSensAllOsc{end,end};
elseif strcmp(modelType,'daAUC')
    outval=1-modelCV.da.ROC{end,end}.AUC;
else
    disp('this is an incorrect modeltype setting')
end
end

