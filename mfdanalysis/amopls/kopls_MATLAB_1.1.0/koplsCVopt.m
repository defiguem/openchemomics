function [modelMain]=koplsCVopt(X,Y,A,oax,modelType,optargin)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for performing K-OPLS cross-validation for a set of
% Y-orthogonal components. The function returns a number of diagnostic
% parameters which can be used to determine the optimal number
% of model components. Optimisation of the kernel parameter is possible
% with grid search within defined limits and simulated annealing.
%
% Use:
% [modelMain]=koplsCVopt(X,Y,A,oax,modelType,{optional,input})
%
% ** INPUT
% X = The measurement matrix.
% Y = The response matrix. Could be binary (for
%   discriminant analysis) or real-valued.
% A = The number of Y-predictive components (integer).
% oax = The number of Y-orthogonal components (integer).
% modelType = 're' for regression, 'da' and 'daAUC' for discriminant 
%   analysis,If 'da' or 'daAUC', sensitivity and specificity will be 
%   calculated together with the area under the ROC curve. For simulated 
%   annealing, optimisation is done for area under ROC curve in 'daAUC' 
%   (only for two-class problems), and for 'da' the mean sensitivity is 
%   optimised.
%
% optargin = cell with optional settings:
%   Possible Input parameters for K-OPLS model, pair-wise input in
%   optargin:
%   kernelType   = kernel type, e.g. 'g' for Gaussian or 'p' for polynomial.
%   preProcK     = Pre-processing settings for the kernel matrix.
%                  Either 'mc' for mean-centering (default) or 'no' for no 
%                  pre-processing.
%   preProcY     = Pre-processing parameter for Y. Either 'mc' for
%                  mean-centering (default), 'uv' for mc + scaling to unit-
%                  variance, 'pa' for mc + Pareto-scaling or 'no' for no 
%                  scaling.
%   opt          = optional settings to do kernel parameter optimisation, 
%                  'no' for no optimisation, uses kernelParams, or 'SA' for
%                  simulated annealing (default for Gaussian kernel), 'GS' 
%                  for gridsearch between values as input in kernelParams 
%                  (default for polynomial kernel)
%   kernelParams = settings for the kernel parameter, leave empty for
%                  optimisation with SA (default), give three values 
%                  [start end nr_steps] for gridsearch, between which the 
%                  gridsearch will be performed and the number of steps; 
%                  note the kernel parameter can not be zero, and polynomial
%                  parameters can only be integer.
%   nrcvinner    = Number of inner cross-validation rounds (integer, 
%                  default = 10).
%   nrcvouter    = Number of outer cross-validation rounds (integer, 
%                  default = 20).
%   cvType       = Type of cross-validation. Either 'nfold' for n-fold
%                  cross-validation, 'mccv' for Monte Carlo CV or 'mccvb' 
%                  for Monte Carlo class-balanced CV; default 'mccv' for 
%                  're' and 'mccvb' for 'da' and 'daAUC' modelType. See 
%                  also 'koplsCrossValSet()' for details.
%   cvFrac       = Fraction of observations (default = 0.75) in the 
%                  training set during cross-validation. Only applicable 
%                  for 'mccv' or 'mccvb' cross-validation (see 'cvType').
%   verbose      = If zero, no output will be displayed (default), if one 
%                  some output will be displayed regarding the 
%                  cross-validation progress. If used with SA optimisation, 
%                  this will also display SA results plot and temperature 
%                  updates.
%   Possible Input parameters for Simulated Annealing, pair-wise input in 
%   optargin:
%   StX       = Starting position for algorithm search - can be initiated
%               automatically or from a number of preset positions. 
%   T_0       = Starting temperature; default = 0.1
%   epsilon   = Termination criterion, default 0.01.
%   Neps      = Number of subsequent optimised points evaluated for 
%               convergence criterion, default 2.
%   v_0       = Starting step vector size, will be adjusted throughout 
%               optimisation but should preferable be able to cover a large 
%               range of the possible optimalkernel parameter values; 
%               default is StX.
%   Ns        = Number of points before reducing vector length, default 5.
%   Nt        = Number of vector reductions before temperature update, 
%               default 5.
%   rT        = Exponential cooling of temperature, default 0.1.
%   nrreps    = Number of runs, default 1.
%
% 
% ** OUTPUT
% modelMain = Object with 'A' predictive components and 'oax'
%   Y-orthogonal components. Contains the following entries:
%  koplsModel = Model calculated with all data and using single (i.e. not-
%     nested) cross-validation optimisation of the kernel parameter.
%  kernelParamslist = (best if optimised) kernel parameter for different
%     CV rounds.
%  KParamfinal = Final kernel parameter used, based on full data set.
%  cv = Cross-validation results:
%    Q2Yhat = Total Q-square result for all Y-orthogonal components.
%    Q2YhatVars = Q-square result per Y-variable for all Y-orthogonal
%       components.
%    Yhat = All predicted Y values as a concatenated matrix.
%    Tcv = Predictive score vector T for all cross-validation rounds.
%    cvTrainIndex = Indices for the training set observations during
%       the cross-validation rounds.
%    cvTestIndex = Indices for the test set observations during the
%       cross-validation rounds.
%  da = Cross-validation results specifically for discriminant
%       analysis (DA) cases:
%    predClass = Predicted class list (integer values).
%    trueClass = True class list (integer values).
%    sensSpec = Sensitivity and specificity values per class and
%       Y-orthogonal components (integer values).
%    confusionMatrix = Confusion matrix during cross-validation
%       rounds.
%    nclasses = Number of classes in model.
%    decisionRule = Decision rule used: 'max' or 'fixed'.
%    ROC = structure containing the AUC (area under ROC curve) and the two
%      vectors based on sensitivity (sens) and 1-specificity (spec) used  
%      for calculation, based on the training data.
%  args = Arguments to the function:
%    A = Number of Y-predictive components.
%    oax = Number of Y-orthogonal components.
%
%  Additional SA output:
%  SAsettings = Settings for simulated annealing.
%
%  Additional GS output:
%  GSresults = optimised predictive performance values for the gridsearched
%  values for the different cross-validation rounds.
%  GSsettings = Grid that was searched.
%  
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
 
release='';
 
[N,m]=size(Y);
 
 
%% % some minor checks....
if isempty(find(strcmp(optargin,'opt'))) 
    if isempty(find(strcmp(optargin,'kernelType'))) 
        kernelType='g';
        opt='SA';   
    else
        kernelType=optargin{(find(strcmp(optargin,'kernelType'))+1)};
        if strcmp(kernelType,'p')
            opt='GS';
            disp('you have to enter three settings (start, end, interval) for gridsearch of a polynomial function, if you have not already done so')
        elseif strcmp(kernelType,'g')
            opt='SA';
        end
    end
else 
    opt=optargin{(find(strcmp(optargin,'opt'))+1)};
    if isempty(find(strcmp(optargin,'kernelType'))) 
        kernelType='g';
    else
        kernelType=optargin{(find(strcmp(optargin,'kernelType'))+1)};
    end
end

if (strcmp(opt,'SA')&strcmp(kernelType,'p'))
    error('use a gridsearch to optimise the polynomial kernel')
end
  

if ~((strcmp(opt,'GS')) | (strcmp(opt,'SA')) | (strcmp(opt,'no')))
    disp('this is an incorrect setting for optimisation, choose ''SA'',''GS'',''no''')
end

if ~((strcmp(kernelType,'g')) | (strcmp(kernelType,'p')))
    disp('this is an incorrect setting for kernelType')
end

if opt~='no'
    tic;
end

if isempty(find(strcmp(optargin,'kernelParams'))) 
    kernelParams=[];
else 
    kernelParams=optargin{(find(strcmp(optargin,'kernelParams'))+1)};
end

if isempty(kernelParams) & strcmp(opt,'no')
    error('if you do not want optimisation of the kernel parameter, enter a value for kernelParams')
elseif strcmp(opt,'GS') & numel(kernelParams)~=3
    error('for gridsearch, enter three values in between which the kernelParameter setting will be chosen, make sure that you only have integer values for polynomial kernels')
end
 
if isempty(find(strcmp(optargin,'nrcvinner')))    
    nrcvinner=10;
else
    nrcvinner=optargin{(find(strcmp(optargin,'nrcvinner'))+1)};
end

if isempty(find(strcmp(optargin,'nrcvouter')))
    nrcvouter=20;
else
    nrcvouter=optargin{(find(strcmp(optargin,'nrcvouter'))+1)};
end

if isempty(find(strcmp(optargin,'cvType')))    
    if strcmp(modelType,'daAUC')|(strcmp(modelType,'da')) 
        cvType='mccvb';
    elseif strcmp(modelType,'re')
        cvType='mccv';
    end
else
    cvType=optargin{(find(strcmp(optargin,'cvType'))+1)};
end

if isempty(find(strcmp(optargin,'cvFrac')))    
    cvFrac=0.75;
else
    cvFrac=optargin{(find(strcmp(optargin,'cvFrac'))+1)};
end

if ((strcmp(modelType,'daAUC'))+(strcmp(modelType,'da'))>0)
   drRule='max'; 
    
   tmp=unique(Y);
   if(all(tmp==[0 1]'))
       if(m==1)
           Y=koplsDummy(Y);
       end
       classVect=koplsReDummy(Y);
       
   elseif(all(mod(Y,1)==0) && m==1)
       classVect=Y;
       Y=koplsDummy(Y+1);
   else
       error('modelType is da, but Y appears to be neither dummy (1 0) matrix nor a vector of (integer) class labels');
   end
   nclasses=length(unique(classVect));
end

if strcmp(cvType,'mccvb') && ~((strcmp(modelType,'da')|strcmp(modelType,'daAUC')))
    error('Class balanced monte-carlo cross validation only applicable to ''da'' or ''daAUC'' modelling');
end
 
if(~any([strcmp(cvType,'mccvb') , strcmp(cvType,'mccv') , strcmp(cvType,'nfold')]))
    error([cvType, '- unknown Cross-validation type']);
end
 
%% convert Y-scaling to more explicit format
if isempty(find(strcmp(optargin,'preProcK'))) 
    preProcK='mc';
else 
    preProcK=optargin{(find(strcmp(optargin,'preProcK'))+1)};
end

if ~((strcmp(preProcK,'mc')) | (strcmp(preProcK,'no')))
    disp('this is an incorrect setting for preProcK scaling')
end

if isempty(find(strcmp(optargin,'preProcY')))
    preProcY='mc';
else 
    preProcY=optargin{(find(strcmp(optargin,'preProcY'))+1)};
end
if ~((strcmp(preProcY,'uv')) | (strcmp(preProcY,'mc')) | (strcmp(preProcY,'no')) | (strcmp(preProcY,'pa')))
    disp('this is an incorrect setting for preProcY scaling')
end

YcenterType='no';
YscaleType='no';
if (~strcmp(preProcY, 'no'))
    YcenterType='mc';
    if (~strcmp(preProcY, 'mc'))
        YscaleType=preProcY;
    end
end

if strcmp(opt,'GS')
    tempKernelParam=kernelParams(1):(kernelParams(2)-kernelParams(1))/(kernelParams(3)-1):kernelParams(2); %change the 10 if you want a different gridsearch sampling scheme
  if strcmp(kernelType,'p') 
    tempKernelParam2=round(tempKernelParam);
    if ~(sum(~(tempKernelParam2==tempKernelParam))==0)
        disp('rounded of your grid search values, please only use integer values with polynomial kernel')
    end
  end
end

if isempty(find(strcmp(optargin,'verbose')))
    verbose=0;
else
    verbose=optargin{(find(strcmp(optargin,'verbose'))+1)};
end


Yhat=ones(size(Y))*NaN;
YhatDaSave=cell(1);
pressyVars=cell(1,1);
pressyVarsTot=cell(1,1);
cvTestIndex=[];
cvTrainingIndex=[];
cvTestIndexbig=[];
cvTrainingIndexbig=[];
 
if (verbose)
    h = waitbar(0,['Please wait... cv round:',num2str(1),' of ',num2str(nrcvouter)]);
end
 
%%
optsetcvvec=[];

for( icv = 1:nrcvouter)

    if (verbose)
        disp(['Cross-validation round: ',num2str(icv),'...']);
        waitbar((icv-1)/nrcvouter,h,['Please wait... cv round:',num2str(icv),' of ',num2str(nrcvouter)]);
    end
 
    %set up CV -----------------
    cvSet=koplsCrossValSet(X,Y,cvFrac,cvType,nrcvouter,icv);
    cvTestIndex=cvSet.testIndex;
    [cvTestIndexbig]=[cvTestIndexbig;cvSet.testIndex];
    cvTrainingIndex=cvSet.trainingIndex;
    cvTrainingIndexbig=[cvTrainingIndexbig;cvSet.trainingIndex];
  
    %----------- opt/SA START -----------------------
    if strcmp(opt,'SA')
        if ~isempty(kernelParams)
            optargin={optargin{1,:},'stX',kernelParams}; 
        end
        [Optset,settings]=koplsSA(X(cvTrainingIndex,:),Y(cvTrainingIndex,:),A,oax,modelType,optargin); 
        modelMain.SAsettings{icv}=settings;    
        modelMain.kernelParamslist{icv}=median(Optset); %median if we did multiple SAs
        kernelParamsfinal=median(Optset); 
    elseif strcmp(opt,'GS')
        for counter=1:length(tempKernelParam)
            modelMain.GSresults{icv}(counter)=koplsModelInternal(X(cvTrainingIndex,:),Y(cvTrainingIndex,:),A,oax,nrcvinner,cvType,preProcK,preProcY,cvFrac,modelType,kernelType,tempKernelParam(counter)); 
        end             
        modelMain.GSsettings{icv}=tempKernelParam;
        [optval,optloc]=min(modelMain.GSresults{icv});
        kernelParamtemp=modelMain.GSsettings{icv}(optloc);
        kernelParamsfinal=kernelParamtemp;
        modelMain.kernelParamslist(icv)=kernelParamtemp;
    elseif strcmp(opt,'no')   
        kernelParamsfinal=kernelParams;
        modelMain.kernelParamslist{icv}=kernelParams;
    else
        error('this is not a valid optimisation setting')
    end
 
               
 
    %----------- SA END -----------------------    
    if kernelParamsfinal==0
        error('note a kernel parameter can not be zero!')
    end
    [K]=koplsKernel(X,[],kernelType,kernelParamsfinal); 
    
    % update the cvSet object with new Kernel parts from optimization
    cvSet.KTrTr=K(cvSet.trainingIndex,cvSet.trainingIndex);
    cvSet.KTeTr=K(cvSet.testIndex,cvSet.trainingIndex);
    cvSet.KTeTe=K(cvSet.testIndex,cvSet.testIndex);
    % ---------------------------------------------------------------          
 
    
    %get Kernel matrices ------- change so that this is done in the K
    %matrix only once and selected by indeces.
    KtrTr=cvSet.KTrTr;
    KteTe=cvSet.KTeTe;
    KteTr=cvSet.KTeTr;
    
    %center Y and kernel matrices ----
    [YScaleObj]=koplsScale(cvSet.yTraining,YcenterType,YscaleType);
    [YScaleObjTest]=koplsScaleApply(cvSet.yTest,YScaleObj);
 
    if (strcmp(preProcK, 'mc'))
        KteTe=koplsCenterKTeTe(KteTe,KteTr,KtrTr);
        KteTr=koplsCenterKTeTr(KteTr,KtrTr);
        KtrTr=koplsCenterKTrTr(KtrTr);
    end
    
    %estimate K-OPLS model-------
    [model]=koplsModel(KtrTr,YScaleObj.X,A,oax,preProcK,preProcY);
    
    %set up model stats----------
    ssy=sum(sum((YScaleObjTest.X).^2));
    ssyVars=(sum((YScaleObjTest.X).^2));    
    ssx=sum(diag(KteTe));
    
    if(icv==1)
        ssyTot=ssy;
        ssyVarsTot=ssyVars;
        ssxTot=ssx;
    else
        ssyTot=ssyTot+ssy;
        ssyVarsTot=ssyVarsTot+ssyVars;        
        ssxTot=ssxTot+ssx;
    end
 
    
    %for each combination of Y-osc components
    for( ioax = 1:oax+1)
        for( ioay = 1:1)       
            [modelPredy]=koplsPredict(KteTr,KteTe,KtrTr, model,ioax-1,0);
            pressy(ioax,ioay)=sum(sum((YScaleObjTest.X-modelPredy.Yhat).^2));        
            pressyVars{ioax,ioay}=(sum((YScaleObjTest.X-modelPredy.Yhat).^2));              
            if((icv==1))
                pressyTot(ioax,ioay)=pressy(ioax,ioay);
                pressyVarsTot{ioax,ioay}=pressyVars{ioax,ioay};
            
            else
                pressyTot(ioax,ioay)=pressyTot(ioax,ioay)+pressy(ioax,ioay);
                pressyVarsTot{ioax,ioay}=pressyVarsTot{ioax,ioay}+pressyVars{ioax,ioay};
             
            end
            
            %if 'da' save Yhat for all rounds
            if ((strcmp(modelType,'daAUC'))+(strcmp(modelType,'da'))>0)
                    if(icv==1)
                        YhatDaSave{ioax,ioay}=[];                        
                    end                    
                    tmp=koplsRescale(YScaleObj,modelPredy.Yhat);
                    YhatDaSave{ioax,ioay}=[YhatDaSave{ioax,ioay};tmp.X];                               
            end
            
            %if highest number of oscs - save Yhat and Xhat
            if(ioax==oax+1)               
                    if(icv==1)
                        Yhat=[];
                    end
                       tmp=koplsRescale(YScaleObj,modelPredy.Yhat);
                       Yhat=[Yhat;tmp.X];
            end
        end
    end
 
    
end %end icv
 
if (verbose)
    waitbar(icv/(nrcvouter),h,'finishing up...');
end
% make final model based on all data 
if opt=='no' 
    KtrTr=K;
    modelMain.koplsModel=koplsModel(KtrTr,Y,A,oax,preProcK,preProcY); 
elseif opt=='SA'
    if (verbose)
    disp('Final model building')
    end
    [Optset,settings]=koplsSA(X,Y,A,oax,modelType,optargin); 
    modelMain.KParamfinal=Optset;
    K=koplsKernel(X,[],kernelType,modelMain.KParamfinal);
    modelMain.koplsModel=koplsModel(K,Y,A,oax,preProcK,preProcY);
    modelMain.time=toc;
elseif opt=='GS'
    if (verbose)
    disp('Final model building')
    end
    for counter=1:length(tempKernelParam)
        ALLGSresults(counter)=koplsModelInternal(X,Y,A,oax,nrcvinner,cvType,preProcK,preProcY,cvFrac,modelType,kernelType,tempKernelParam(counter)); 
    end
    ALLGSsettings=tempKernelParam;
    [ALLoptval,ALLoptloc]=min(ALLGSresults);
    modelMain.KParamfinal=ALLGSsettings(ALLoptloc);
    K=koplsKernel(X,[],kernelType,modelMain.KParamfinal);
    modelMain.koplsModel=koplsModel(K,Y,A,oax,preProcK,preProcY);
    modelMain.time=toc;
end
    
modelMain.cv.Yhat=Yhat;
modelMain.cv.Tcv=Yhat*modelMain.koplsModel.Cp*modelMain.koplsModel.Bt{oax+1};
modelMain.cv.Q2Yhat=[]; 
modelMain.cv.Q2YhatVars=cell(1,1);
 
    for( ioax = 1:oax+1)
        for( ioay = 1:1) 
            modelMain.cv.Q2Yhat(ioax,ioay)=1-pressyTot(ioax,ioay)./ssyTot;
            modelMain.cv.Q2YhatVars{ioax,ioay}=1-pressyVarsTot{ioax,ioay}./ssyVarsTot;
        end
    end
 
    modelMain.cv.cvTestIndex=cvTestIndexbig;
    modelMain.cv.cvTrainingIndex=cvTrainingIndexbig;
 
    if ((strcmp(modelType,'daAUC'))+(strcmp(modelType,'da'))>0)
 
        %get sens/spec for each y-orth component... eval of model
        for( i = 1:oax+1) %we would have no osc comps for dummy matrix...
                if(strcmp(drRule,'max'))
                    predClass=koplsMaxClassify(YhatDaSave{i,1});
                elseif(strcmp(drRule,'fixed'))
                    predClass=koplsBasicClassify(YhatDaSave{i,1},1/nclasses);
                else
                    warning(['Decision rule given: ',drRule,' is not valid/implemnted'])
                end
                [da.sensAllOsc{i}, da.specAllOsc{i}, da.classvecAllOsc{i}, da.tot_sensAllOsc{i},da.meanSensAllOsc{i},da.meanSpecAllOsc{i}]=koplsSensSpec(classVect(cvTestIndexbig), predClass);
                  if (size(YhatDaSave{i,1},2))<3&&(length(unique(classVect(cvTestIndexbig)))==2) 
                    [da.ROC{i}]=koplsRoc(YhatDaSave{i,1},classVect(cvTestIndexbig));
                end
        end

        % get sens/spec for max number of oscs
        if(strcmp(drRule,'max'))
            predClass=koplsMaxClassify(Yhat);
        elseif(strcmp(drRule,'fixed'))
            predClass=koplsBasicClassify(Yhat,1/nclasses);
        else
            warning(['Decision rule given: ',drRule,' is not valid/implemnted'])
        end
        
           [da.sens, da.spec, da.classvec, da.tot_sens,da.meanSens,da.meanSpec]=koplsSensSpec(classVect(cvTestIndexbig), predClass);
           [da.confusionMatrix]=koplsConfusionMatrix(classVect(cvTestIndexbig), predClass);
           da.trueClass=classVect(cvTestIndexbig);
           da.nclasses=nclasses;
        modelMain.da=da;
        modelMain.da.predClass=predClass;        
        modelMain.da.decisionRule=drRule;
                %CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
                %compatibility and comparison w/ simca-p etc
        if(strcmp(cvType,'nfold'))
            [tmp,cvOrder]=sort(cvTestIndex);
            modelMain.da.predClass=modelMain.da.predClass(cvOrder);
            modelMain.da.trueClass=modelMain.da.trueClass(cvOrder);           
        end
        
    end
 
    %CHANGE TO ORIGINAL ORDER IF NFOLD CV - for backward
    %compatibility and comparison w/ simca-p etc
    if(strcmp(cvType,'nfold'))
        [tmp,cvOrder]=sort(cvTestIndexbig);
            modelMain.cv.Yhat=modelMain.cv.Yhat(cvOrder,:);
           modelMain.cv.Tcv=modelMain.cv.Tcv(cvOrder,:);
 
    end
 
if (verbose)
    close(h);
end
 
modelMain.release=release;
modelMain.args.oax=oax;
modelMain.args.A=A;
modelMain.class='koplscvopt';
