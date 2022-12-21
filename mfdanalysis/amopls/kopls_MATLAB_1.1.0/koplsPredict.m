function [modelp]=koplsPredict(KteTr,Ktest,Ktrain,model,nox,rescaleY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Performs prediction of new samples from an existing K-OPLS model
% (see koplsModel()' to calculate K-OPLS models).
% The function projects the Y-predictive and Y-orthogonal scores
% components to predict a value of the response matrix Y.
% The dimensionality of the parameters is determined from the
% specified model.
%
% ** INPUT
% KteTr = The hybrid test/training kernel matrix;
%   KteTr = <phi(Xte),phi(Xtr)>.
% Ktest = The pure test kernel matrix;
%   Ktest = <phi(Xte),phi(Xte)>.
% Ktrain = The training kernel matrix (same as used in
%   model training); Ktrain = <phi(Xtr),phi(Xtr)>.
% model = K-OPLS model object.
% nox = Number of Y-orthogonal components. If not specified, the
%   number used during model training will be employed.
% rescaleY = Boolean parameter. If true, predicted values of the
%   response (Yhat) is rescaled according to the pre-processing
%   settings of the model. If false, Yhat is not rescaled (default).
%
% ** OUTPUT:
% modelp = Object with the following entries:
%   Tp = Predicted predictive score matrix for all generations
%       0:'nox' of Y-orthogonal vectors.
%   T = Predictive score matrix for the final model with 'nox'
%       Y-orthogonal vectors.
%   to = Predicted Y-orthogonal score vectors.
%   EEprime = Calculated residuals for the test kernel 'Ktest',
%       useful e.g. for residual statistics.
%   Yhat = Predicted values of the response matrix.
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


if (~strcmp(model.class,'kopls'))
    error('Model must be of type "kopls". Aborting.')
end


%% mean centering of K matrices

%% centering (order of these arg is important...)
KteTeMc=Ktest;
if (strcmp(model.preProc.K,'mc'))
	KteTeMc=koplsCenterKTeTe(Ktest,KteTr,Ktrain);
end
KteTe=cell(model.nox+1,model.nox+1);
KteTe{1,1}=KteTeMc;

KteTrMc=KteTr;
if (strcmp(model.preProc.K,'mc'))
	KteTrMc=koplsCenterKTeTr(KteTr,Ktrain);
end

KteTr=cell(model.nox+1,model.nox+1);
KteTr{1,1}=KteTrMc;

%% init of Y-orth comps
to=[];

%% check if last arg is number of components to use in prediction:

if(nargin>=5)
    if(nox>model.nox)
        warning('Number of Y-orthogonal components to use is higher than in model - setting number of Yorth to max in model');
        nox=model.nox;
    end
else
    nox=model.nox;
end

if (nargin < 6)
    rescaleY=0;
end

%% KOPLS prediction


for(i=1:nox) %step1
    
    %step2
    Tp{i}=KteTr{i,1}*model.Up*model.Sp^(-1/2);
    %Yhat{i}=Tp{i}*model.Bt{i}*model.Cp';

    %step3
    to{i}=(KteTr{i,i}-Tp{i}*model.Tp{i}')*model.Tp{i}*model.co{i}*model.so{i}^(-1/2);
    
    %step4
    to{i}=to{i}./model.toNorm{i};

    % step 4.5 deflate KteTe. (this is an EXTRA feature - not in alg. in
    % paper )
    KteTe{i+1,i+1} = KteTe{i,i} - KteTr{i,i}*model.to{i}*to{i}' - to{i}*model.to{i}'*KteTr{i,i}' + to{i}*model.to{i}'*model.K{i,i}*model.to{i}*to{i}';
    
    %step5
    KteTr{i+1,1}=KteTr{i,1} -to{i}*model.to{i}'*model.K{1,i}';
    
    %step6
    KteTr{i+1,i+1}=KteTr{i,i}-KteTr{i,i}*model.to{i}*model.to{i}'-to{i}*model.to{i}'*model.K{i,i}+to{i}*model.to{i}'*model.K{i,i}*model.to{i}*model.to{i}';
    
end 
    %step7

 if(nox==0)
     i=0;
 end
    
Tp{i+1}=KteTr{i+1,1}*model.Up*model.Sp^(-1/2);
%Yhat{i+1}=Tp{i+1}*model.Bt{i+1}*model.Cp';
Yhat=Tp{i+1}*model.Bt{i+1}*model.Cp';

if (rescaleY)
	
    if (~strcmp(model.preProc.Y,'no'))
        YhatRescaled=koplsRescale(model.preProc.paramsY, Yhat);
        Yhat=YhatRescaled.X;
    else
        warning('Attempted re-scale of Yhat although no pre-processing parameters have been set.')
    end
end



%---- Extra stuff ----------------------------------
%this appears to be correct - but does not match previous code...
EEprime=(KteTe{i+1,i+1}-Tp{i+1}*Tp{i+1}'); 
%--------------------------------------------------

modelp.Tp=Tp;
modelp.to=to;
%modelp.KteTr=KteTr;
modelp.EEprime=EEprime;
modelp.Yhat=Yhat;

end

