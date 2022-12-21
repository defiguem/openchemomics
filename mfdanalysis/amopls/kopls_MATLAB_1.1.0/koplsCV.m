function [modelMain]=koplsCV(K,Y,A,oax,nrcv,cvType,preProcK,preProcY,cvFrac,modelType,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for performing K-OPLS cross-validation for a set of
% Y-orthogonal components. The function returns a number of diagnostic
% parameters which can be used to determine the optimal number
% of model components.
%
% ** INPUT
% K = The kernel matrix (un-centered); see 'koplsKernel()' for details.
% Y = The response matrix (un-centered/scaled). Could be binary (for
%   discriminant analysis) or real-valued.
% A = The number of Y-predictive components (integer).
% oax = The number of Y-orthogonal components (integer).
% nrcv = Number of cross-validation rounds (integer).
% cvType = Type of cross-validation. Either 'nfold' for n-fold
%   cross-validation, 'mccv' for Monte Carlo CV or 'mccvb' for
%   Monte Carlo class-balanced CV. See also 'koplsCrossValSet()' for
%   details.
% preProcK = Pre-processing settings for the kernel matrix.
%   Either 'mc' for mean-centering or 'no' for no pre-processing.
% preProcY = Pre-processing parameter for Y. Either 'mc' for
%   mean-centering, 'uv' for mc + scaling to unit-variance,
%   'pa' for mc + Pareto-scaling or 'no' for no scaling.
% cvFrac = Fraction of observations in the training set during
%   cross-validation. Only applicable for 'mccv' or 'mccvb'
%   cross-validation (see 'cvType').
% modelType = 'da' for discriminant analysis, 're' for regression.
%   If 'da', sensitivity and specificity will be calculated.
% verbose = If zero, no output will be displayed, otherwise some
%   output will be displayed regarding the cross-validation progress
% 	(default).
%
% ** OUTPUT
% modelMain = Object with 'A' predictive components and 'oax'
%   Y-orthogonal components. Contains the following entries:
%  cv = Cross-validation results:
%  	 Q2Yhat = Total Q-square result for all Y-orthogonal components.
%	 Q2YhatVars = Q-square result per Y-variable for all Y-orthogonal
%       components.
%	 Yhat = All predicted Y values as a concatenated matrix.
%	 Tcv = Predictive score vector T for all cross-validation rounds.
%	 cvTrainIndex = Indices for the training set observations during
%       the cross-validation rounds.
%	 cvTestIndex = Indices for the test set observations during the
%       cross-validation rounds.
%  da = Cross-validation results specifically for discriminant
%       analysis (DA) cases:
%    predClass = Predicted class list per class and Y-orthogonal
%       components (integer values).
%	 trueClass = Predicted class list per class and Y-orthogonal
%       components (integer values).
%	 sensSpec = Sensitivity and specificity values per class and
%       Y-orthogonal components (integer values).
%	 confusionMatrix = Confusion matrix during cross-validation
%       rounds.
%	 nclasses = Number of classes in model.
%	 decisionRule = Decision rule used: 'max' or 'fixed'.
%  args = Arguments to the function:
%  	 A = Number of Y-predictive components.
%	 oax = Number of Y-orthogonal components.
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

if (nargin < 11)
    verbose=1;
end


%%TODO: set default values for the last entries (nrcv/cvType/etc.)


%%
% some minor checks....
    if ((strcmp(modelType,'daAUC'))+(strcmp(modelType,'da'))>0)
   drRule='max'; %move to arg... %this is a prameter for DA decision rule
    
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

if(strcmp(cvType,'mccvb') && ~strcmp(modelType,'da'))
    error('Class balanced monte-carlo cross validation only applicable to da modelling');
end

if(~any([strcmp(cvType,'mccvb') , strcmp(cvType,'mccv') , strcmp(cvType,'nfold')]))
    error([cvType, '- unknown Cross-validation type']);
end
%%


%% convert Y-scaling to more explicit format
YcenterType='no';
YscaleType='no';
if (~strcmp(preProcY, 'no'))
    YcenterType='mc';
    if (~strcmp(preProcY, 'mc'))
        YscaleType=preProcY;
    end
end

%keyboard;


Yhat=ones(size(Y))*NaN;
YhatDaSave=cell(1);
pressyVars=cell(1,1);
pressyVarsTot=cell(1,1);
cvTestIndex=[];
cvTrainingIndex=[];

if (verbose)
    h = waitbar(0,['Please wait... cv round:',num2str(1),' of ',num2str(nrcv)]);
end

%%
for( icv = 1:nrcv)
    
    if (verbose)
        disp(['Cross-validation round: ',num2str(icv),'...']);
        waitbar((icv-1)/nrcv,h,['Please wait... cv round:',num2str(icv),' of ',num2str(nrcv)]);
    end

    %set up CV -----------------
	cvSet=koplsCrossValSet(K,Y,cvFrac,cvType,nrcv,icv);
    cvTestIndex=[cvTestIndex;cvSet.testIndex];
    cvTrainingIndex=[cvTrainingIndex;cvSet.trainingIndex];
    
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
    [model]=koplsModel(KtrTr,YScaleObj.X,A,oax,'no','no');
    
    %set up model stats----------
    ssy=sum(sum((YScaleObjTest.X).^2));
    ssyVars=(sum((YScaleObjTest.X).^2));	
    ssx=sum(diag(KteTe));
    
    if(icv==1)
        ssyTot=ssy;
        ssyVarsTot=ssyVars;
        ssxTot=ssx;
        %ssxVarsTot=ssxVars;
    else
        ssyTot=ssyTot+ssy;
        ssyVarsTot=ssyVarsTot+ssyVars;        
        ssxTot=ssxTot+ssx;
        %ssxVarsTot=ssxVarsTot+ssxVars;
    end
    
    
    %for each combination of Y-osc components
	for( ioax = 1:oax+1)
        for( ioay = 1:1)%oay+1)        
            % KOPLS predict yYhat                                   
            [modelPredy]=koplsPredict(KteTr,KteTe,KtrTr, model,ioax-1,0);
           
           

            pressy(ioax,ioay)=sum(sum((YScaleObjTest.X-modelPredy.Yhat).^2));        
            pressyVars{ioax,ioay}=(sum((YScaleObjTest.X-modelPredy.Yhat).^2));	            
            
            %keyboard;
            
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
                 
                    
                    %+mean on Yhat
                    tmp=koplsRescale(YScaleObj,modelPredy.Yhat);
                    YhatDaSave{ioax,ioay}=[YhatDaSave{ioax,ioay};tmp.X];                               
            end
            
            %if highest number of oscs - save Yhat and Xhat
            if(ioax==oax+1)% && ioay==oay+1)                
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
    waitbar(icv/(nrcv),h,'finishing up...');
end

%[scaleY]=koplsScale(Y,YcenterType,YscaleType);
%KtrTr=koplsKernel(X,X,kernelType,kernelParam);
%if (strcmp(preProcK,'mc'))
%    KtrTr=koplsCenterKTrTr(KtrTr);
%end
%modelMain.koplsModel=koplsModel(KtrTr,scaleY.X,A,oax,preProcK,preProcY);

KtrTr=K;
modelMain.koplsModel=koplsModel(KtrTr,Y,A,oax,preProcK,preProcY);

modelMain.cv.Yhat=Yhat;

%is this correct?
modelMain.cv.Tcv=Yhat*modelMain.koplsModel.Cp*modelMain.koplsModel.Bt{oax+1};

    modelMain.cv.Q2Yhat=[]; 
    modelMain.cv.Q2YhatVars=cell(1,1);


    for( ioax = 1:oax+1)
        for( ioay = 1:1)%oay+1)   
            modelMain.cv.Q2Yhat(ioax,ioay)=1-pressyTot(ioax,ioay)./ssyTot;
            modelMain.cv.Q2YhatVars{ioax,ioay}=1-pressyVarsTot{ioax,ioay}./ssyVarsTot;
        end
    end

    modelMain.cv.cvTestIndex=cvTestIndex;
    modelMain.cv.cvTrainingIndex=cvTrainingIndex;

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
                %keyboard;
                [da.sensAllOsc{i}, da.specAllOsc{i}, da.classvecAllOsc{i}, da.tot_sensAllOsc{i},da.meanSensAllOsc{i},da.meanSpecAllOsc{i}]=koplsSensSpec(classVect(cvTestIndex), predClass);
                if (size(YhatDaSave{i,1},2))<3&&(length(unique(classVect(cvTestIndex)))==2) 
                [da.ROC{i}]=koplsRoc(YhatDaSave{i,1},classVect(cvTestIndex));
                end
        end
        

        
        % get sens/spec for max number of oscs.... (hmm redundant).
        
        if(strcmp(drRule,'max'))
            predClass=koplsMaxClassify(Yhat);
        elseif(strcmp(drRule,'fixed'))
            predClass=koplsBasicClassify(Yhat,1/nclasses);
        else
            warning(['Decision rule given: ',drRule,' is not valid/implemnted'])
        end
        

           [da.sens, da.spec, da.classvec, da.tot_sens,da.meanSens,da.meanSpec]=koplsSensSpec(classVect(cvTestIndex), predClass);
           [da.confusionMatrix]=koplsConfusionMatrix(classVect(cvTestIndex), predClass);
           da.trueClass=classVect(cvTestIndex);
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


    %CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
    %compatibility and comparison w/ simca-p etc
    if(strcmp(cvType,'nfold'))
        [tmp,cvOrder]=sort(cvTestIndex);

           modelMain.cv.Yhat=modelMain.cv.Yhat(cvOrder,:);
           modelMain.cv.Tcv=modelMain.cv.Tcv(cvOrder,:);

    end

if (verbose)
    close(h);
end

modelMain.release=release;
modelMain.args.oax=oax;
%modelMain.args.oay=oay;
modelMain.args.A=A;
modelMain.class='koplscv';

end






