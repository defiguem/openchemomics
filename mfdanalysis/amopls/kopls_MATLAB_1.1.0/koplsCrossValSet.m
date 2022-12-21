function [cvSet]=koplsCrossValSet(K,Y,modelFrac,type,nfold,nfoldRound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates sets of training/test observations useful for cross-
% validation (CV). How the sets are generated is determined by the
% 'type' parameter, which can be either 'nfold' for n-fold cross-
% validation, 'mccv' for Monte Carlo CV, 'mccvb' for Monte Carlo
% class-balanced CV. 
%
% ** INPUT
% K = Kernel matrix.
% Y = Response matrix.
% type = Type of cross-validation:
%	'nfold' for n-fold, 'mccv' for Monte Carlo CV, 'mccvb' for
%	Monte Carlo class-balanced CV.
%  nfold = Number of total nfold rounds (if type='nfold').
%  nfoldRound = Current nfold round (if type='nfold').
%
% ** OUTPUT
%  Object 'cvSet' with the following entries:
%   type = Cross-validation type.
%   nfold = Number of nfold rounds.
%   nfoldRound = The current nfold round.
%   KTrTr = Kernel training matrix; KTrTr = <phi(Xtr),phi(Xtr)>.
%   KTeTr = Kernel test/training matrix;
%       KTeTr = <phi(Xte),phi(Xtr)>.}
%   KTeTe = Kernel test matrix; KTeTe = <phi(Xte),phi(Xte)>.
%   yTraining = Y training set.
%   yTest = Y test set.
%   trainingIndex = Indices of training set observations.
%   testIndex = Indices of test set observations.
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


    modInd=[];
    predInd=[];

	
    if(strcmp(type,'mccvb'))  %'Monte-carlo cv - class Balanced'
            %check if Y is dummy or labels...
            tmp=unique(Y);
            if(all(tmp==[0 1]'))
                classVect=koplsReDummy(Y);
            else
                classVect=Y;
            end
                
		minset=unique(classVect); %find all classlabels
        for i=1:length(minset) %for each class
            currentClass=minset(i); %current class label
            ind=find(classVect==currentClass); %find all samples from current class

            %randomize
            ran=rand(length(ind),1); %randomize
            [tmp,rand_ind]=sort(ran); %sort randomize number to get randomized index string
            ind=ind(rand_ind); %apply randomization on the real index vector
            %-end randomize

            modelLim=ceil(length(ind)*modelFrac); %number of elemnts that will get assigned as model
            modInd=[modInd;ind(1:modelLim)];
            predInd=[predInd;ind(modelLim+1:end)];
        end
    end
	
	
    if(strcmp(type,'mccv'))  %'Monte-carlo cv'		
            %randomize
            ran=rand(length(K(:,1)),1); %randomize
            [tmp,rand_ind]=sort(ran); %sort randomize number to get randomized index string
            ind=[1:length(ran)]';
            ind=ind(rand_ind); %apply randomization on the real index vector
            modelLim=ceil(length(ind)*modelFrac); %number of elemnts that will get assigned as model
            modInd=[ind(1:modelLim)];
            predInd=[ind(modelLim+1:end)];        
    end	
	
    if(strcmp(type,'nfold'))  %'N-Fold cross validation'
        predInd=[nfoldRound:nfold:length(Y(:,1))]';
        modInd=[setdiff(1:length(Y(:,1)),predInd)]';
    end

	[a,b]=size(K);

    cvSet.type=type;
    cvSet.nfold=nfold;
    cvSet.nfoldRound=nfoldRound;

	if(a==b)
		cvSet.KTrTr=K(modInd,modInd);
		cvSet.KTeTr=K(predInd,modInd);
		cvSet.KTeTe=K(predInd,predInd);
	else	
		cvSet.KTrTr=[];
		cvSet.KTeTr=[];
		cvSet.KTeTe=[];
	end
	 
    cvSet.yTraining=Y(modInd,:);
    cvSet.yTest=Y(predInd,:);

    cvSet.trainingIndex=modInd;
    cvSet.testIndex=predInd;
        
end
