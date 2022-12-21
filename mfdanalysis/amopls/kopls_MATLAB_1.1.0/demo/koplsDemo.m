%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script contains a demonstration of the functionality in the
% 'kopls' package using a simulated data set. The data set is
% represented by 1000 spectral variables from two different classes
% and is available in the 'koplsDemo.mat' workspace. The
% demonstration essentially consists of two main steps.
%
% The first step is to demonstrate how K-OPLS handles the
% model evaluation (using cross-validation), model building and
% subsequent classification of external data from a non-linear data
% set. The second step is to demonstrate how K-OPLS works in the
% presence of response-independent (Y-orthogonal) variation, using
% the same data set but with a strong systematic class-specific
% disturbance added.
%
% ** THE 'koplsDemo.mat' WORKSPACE
% The koplsDemo.mat workspace contains the following objects:
%   Xtr = The training data matrix, with 400 observations and
%       1000 spectral variables.
%   Xte = The test data matrix, with 400 observations and
%       1000 spectral variables.
%   Xtro = Same data as 'Xtr', but with class-specific systematic
%       noise added.
%   Xteo = Same data as 'Xte', but with class-specific systematic
%       noise added.
%   Ytr = A binary matrix of class assignments for the training data.
%   Yte = A binary matrix of class assignments for the test data.
%   class1 = A vector of indicies of the samples in class #1.
%   class2 = A vector of indicies of the samples in class #2.
% 
%
% ** INSTRUCTIONS
% 1) Make sure that all 'kopls*.m' files are in the current path.
% 2) Run 'koplsDemo'.
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

echo on

%% Make sure that the workspace is loaded, otherwise try to load it
try
   Xtr; Xte; Xtro; Xteo; Ytr; Yte; class1; class2;
catch
    disp('Demo workspace not loaded, trying to load it automatically...');
    
    try
        load('koplsDemo.mat');
        Xtr; Xte; Xtro; Xteo; Ytr; Yte; class1; class2;
        disp('Successfully loaded workspace koplsDemo.mat');
    catch
        error('Demo workspace has to be loaded manually before running the demo; please consult the documentation for details');
    end
end


%%%%%%%% START OF DEMO

%% Now running Principal Component Analysis (PCA) on simulated training
%% data
[U,D,V]=svds(Xtr, 2);
  
%% Plot PCA score vectors to demonstrate data set properties, also
%% demonstrate 'time' dependency of data
figure;
plot(U(class1,1), U(class1,2), 'bx');
hold on;
plot(U(class2,1), U(class2,2), 'ro');
title('PCA of original data set');
hold off;

%%%%%%%%%%%%%%%% K-OPLS modeling

%% Define kernel function parameter
%% This sets the kernel parameter for gaussian kernel to 25, at the end of this demo it is shown how to optimise this parameter.
sigma=25;

%% Construct the training kernel
Ktr=koplsKernel(Xtr,[],'g',sigma);

%% Find optimal number of Y-orthogonal components by cross-validation
%% This step performs the CV
modelCV=koplsCV(Ktr,Ytr,1,3,7,'nfold','mc','mc',0.75,'da');

%% Plot cross-validation results
figure;
koplsPlotCVDiagnostics(modelCV);

%% Plot sensitivity and specificity measures
figure;
koplsPlotSensSpec(modelCV);
title('Sensitivity and specificity results during cross-validation');

%% 'nox' defines the number of Y-orthogonal components used in the final
%% model. This value is selected according to the cross-validation
%% results. We pick nox=1 for diagnostics, although the model with
%% nox=0 is more or less equivalent.
nox=1;

%% Create test kernels, using the same function as the test data
%% but with now involving the test set matrix
KteTr=koplsKernel(Xte,Xtr,'g',sigma);
KteTe=koplsKernel(Xte,[],'g',sigma);

%% Construct final model
modelOrg=koplsModel(Ktr,Ytr,1,nox,'mc','mc');

%% Predict test set
modelOrgPred=koplsPredict(KteTr,KteTe,Ktr,modelOrg,nox,1);

%% View scores from the final model
figure;
koplsPlotScores(modelOrg)

%% View predictions for external test set
figure;
plot(modelOrgPred.Yhat, Yte, 'xb');
hold on;
xlabel('Predicted');
ylabel('Observed');
title('Obs. vs. pred. for original data');
plot( [0.5 0.5], [0 1], 'r--');
hold off;

%%%%%%%%%%%%%%%%%%% Now model data with Y-orthogonal variation added (to one class)


%% Now running Principal Component Analysis (PCA) on simulated training
%% data with one class distorted
[Uo,Do,Vo]=svds(Xtro, 2);

%% Plot PCA score vectors to demonstrate data set properties
figure;
plot(Uo(class1,1), Uo(class1,2), 'bx')
hold on;
plot(Uo(class2,1), Uo(class2,2), 'ro')
title('PCA of data with Y-ortho variation added')
hold off

%% Create training and test kernels
Ktro=koplsKernel(Xtro,[],'g',sigma);
KteTro=koplsKernel(Xteo,Xtro,'g',sigma);
KteTeo=koplsKernel(Xteo,[],'g',sigma);

%% Model and predict
modelOSC=koplsModel(Ktro,Ytr,1,nox,'mc','mc');
modelOSCPred=koplsPredict(KteTro,KteTeo,Ktro,modelOSC,nox,1);

%% We skip cross-validation this time (for speed)
figure;
koplsPlotModelDiagnostics(modelOSC);

%% Plot the scores. Note bimodality in predictive component and
%% 'trimodality' in first Y-ortho component.
figure;
koplsPlotScores(modelOSC);

%% View predictions for external test set
figure;
plot(modelOSCPred.Yhat, Yte, 'bx');
hold on
xlabel('Predicted');
ylabel('Observed');
title('Obs. vs. pred. with Y-ortho variation added');
plot( [0.5 0.5], [0 1], 'r--');
hold off;


%% Time analysis to demonstrate SA; as this takes some time we only use a fraction of the data:
Xtr2=Xtr(1:10:400,:);
Ytr2=[1:20,1:20]';
[U2,D2,V2]=svds(Xtr2, 2);
  
%% Plot PCA score vectors to demonstrate data set properties, also
%% demonstrate 'time' dependency of data
figure;
plot(U2(1:20,1), U2(1:20,2), 'bx');hold on; 
hold on;
plot(U2(21:40,1), U2(21:40,2), 'ro'); text(U2(:,1), U2(:,2),num2str(Ytr2))
title('PCA of part original data set with time indicated');
hold off;

%% try and predict the timepoint Ytimetr, optimise the kernel parameter 
%% settings and model using gridsearch between 1 and 30 in 30 steps, 
%% gaussian kernel (default) is used
[modelMain_GSre]=koplsCVopt(Xtr2,Ytr2,1,3,'re',{'kernelParams',[1 30 30],'opt','GS','nrcvouter',5}); 
%% to plot the gridsearch results
koplsPlotOptResults(modelMain_GSre,'GS','re');

%% Optimise the kernel parameter settings and model using simulated
%% annealing, fast cooling and a Gaussian kernel are used with 5 outer
%% cross validation loops (nrcvouter), and 10 inner cross-validation 
%% loops (default), 'verbose',1 will display the cross-validation and
%% SA optimisation schedule.
[modelMain_SAre]=koplsCVopt(Xtr2,Ytr2,1,3,'re',{'verbose',1,'t_0',10,'rT',0.2,'opt','SA','nrcvouter',5,'kernelType','g'}); 

%% to plot the simulated annealing results
koplsPlotOptResults(modelMain_SAre,'SA','re');
% note that not all evaluated points are plotted, for clarity.

%% other test data set
Xtr3=Xtro([1:10:400],:);
Ytr3=[1:20,1:20]';
%% gridsearch with a polynomial function
[modelMain_GSpoly]=koplsCVopt(Xtr3,Ytr3,1,4,'re',{'kernelParams',[1 5 5],'opt','GS','kernelType','p','cvType','nfold'}); 

%% to plot the gridsearch results
koplsPlotOptResults(modelMain_GSpoly,'GS','re');

%% Similar for optimisation of a discriminant analysis problem, using area
%% under ROC curve and based on the distorted the data;
%% simulated annealing with (20 outer loops so might take a few minutes)
Xtr4=Xtro([1:15,201:215],:);
Ytr4=Ytr([1:15,201:215],:);
[modelMain_SAda]=koplsCVopt(Xtr4,Ytr4,1,1,'daAUC',{}); 

%% plot of the optimisation: 
koplsPlotOptResults(modelMain_SAda,'SA','daAUC');

%% a plot of the predictive and orthogonal scores:
figure;hold on
plot(modelMain_SAda.koplsModel.T([1:15]),modelMain_SAda.koplsModel.To([1:15]),'bx')
plot(modelMain_SAda.koplsModel.T([16:30]),modelMain_SAda.koplsModel.To([16:30]),'ro')
xlabel('Predictive component');
ylabel('Orthogonal Component');
title('Score plot of SA-K-OPLS discriminant analysis model');
legend('class 1','class 2','Location','Best')

%%%%%%%% END OF DEMO

echo off

