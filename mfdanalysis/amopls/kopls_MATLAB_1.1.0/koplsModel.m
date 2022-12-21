function [model]=koplsModel(K,Y,A,nox,preProcK,preProcY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for training a K-OPLS model. The function constructs a
% predictive regression model for predicting the values of 'Y by
% using the information in 'K'. The explained variation is separated
% into predictive components (dimensionality is determined by the
% parameter 'A') and 'Y'-orthogonal components (dimensionality
% determined by the parameter 'nox').
%
% ** INPUT
% K = Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
% Y = Response matrix (un-centered/scaled).
% A = Number of predictive components.
% nox = Number of Y-orthogonal components.
% preProcK = Pre-processing parameters for the 'K' matrix:
%   'mc' for mean-centering, 'no' for no centering.
% preProcY = Pre-processing parameters for the 'Y' matrix:
%   'mc' for mean-centering, 'uv' for mc + scaling to unit variance,
%   'pa' for mc + Pareto, 'no' for no scaling.
%
% ** OUTPUT
% model = Object with the following entries:
%   Cp = Y loading matrix.
%   Sp = Sigma matrix, containing singular values from Y'*K*Y used
%       for scaling.
%   Sps = Sp^(-1/2).
%   Up = Y score matrix.
%   Tp = Predictive score matrix for all Y-orthogonal components.
%   T = Predictive score matrix for the final model.
%   co = Y-orthogonal loading vectors.
%   so = Eigenvalues from estimation of Y-orthogonal loading vectors.
%   To = Y-orthogonal score matrix.
%   toNorm = Norm of the Y-orthogonal score matrix prior to scaling.
%   Bt = T-U regression coefficients for predictions.
%   A = Number of predictive components.
%   nox = Number of Y-orthogonal components.
%   K = The kernel matrix.
%   EEprime = The deflated kernel matrix for residual statistics.
%   sstot_K = Total sums of squares in 'K'.
%   R2X = Cumulative explained variation for all model components.
%   R2XO = Cumulative explained variation for Y-orthogonal
%       model components.
%   R2XC = Explained variation for predictive model components after
%       addition of Y-orthogonal model components.
%   sstot_Y = Total sums of squares in Y.
%   R2Y = Explained variation of Y.
%   preProc = Pre-processing parameters: 
%  	  K = Pre-processing setting for K = 'preProcK'.
%	  Y = Pre-processing setting for Y = 'preProcY'.
%	  paramsY = Scaling parameters for Y.
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


[nn,kk]=size(K);
I=eye(nn);
Kmc=K;
if (strcmp(preProcK,'mc'))
    [Kmc]=koplsCenterKTrTr(K);
end
K=cell(nox+1,nox+1);
K{1,1}=Kmc;


Y_old=Y;
scaleParams=cell(1);
if (strcmp(preProcY,'mc') || strcmp(preProcY,'uv') || strcmp(preProcY,'pa'))
	
    if(strcmp(preProcY, 'mc'))
        scale='no';
    else
        scale=preProcY;
    end
    
    scaleParams=koplsScale(Y, 'mc', scale);
	Y=scaleParams.X;
end




%mean centering ---------------------
% I=eye(length(K(:,1)));
% I_n=ones(length(K(:,1)),1);
% n=length(K(:,1));
% 
% %this is the mean center step in feature space..
% Kmc = (I- (1/n).* I_n*I_n') * K*(I-(1/n).*I_n*I_n');
% %Kold_mc=K;




%% initiate Yorth related vars
to=[];
co=[];
so=[];
toNorm=[];

%% KOPLS mode estimation -------------

%step 1
[Cp,Sp,V] = svd(Y'*K{1,1}*Y);
%[Cp,Sp]=eigs(K,A);
Cp=Cp(:,1:A);
Sp=Sp(1:A,1:A);

%step2
Up=Y*Cp;


for(i=1:nox)%Step3

    %step4
    Tp{i}=K{1,i}'*Up*Sp^(-1/2);
    Bt{i}=inv(Tp{i}'*Tp{i})*Tp{i}'*Up;  
    
    %step5
    [CoTmp,SoTmp,Vo] = svd( Tp{i}'* (K{i,i}-Tp{i}*Tp{i}')* Tp{i} );
    %[Co{i},So{i}]=eigs( T{i}'* (K{i,i}-T{i}*T{i}')* T{i} );
    co{i}=CoTmp(:,1);
    so{i}=SoTmp(1,1);
    
    %step6
    to{i}=(K{i,i}-Tp{i}*Tp{i}')*Tp{i}*co{i}*so{i}^(-1/2);
    
    %step7
    toNorm{i}=sqrt(to{i}'*to{i});
    
    %step8
    to{i}=to{i}./toNorm{i};
    
    %step9
    K{1,i+1}=K{1,i}*(I - to{i}*to{i}');
    
    %step10
    K{i+1,i+1}=(I - to{i}*to{i}')*K{i,i}*(I - to{i}*to{i}');
    
    
     
end %step 11

%step12
Tp{nox+1}=K{1,nox+1}'*Up*Sp^(-1/2);

%Step13
Bt{i+1}=inv(Tp{nox+1}'*Tp{nox+1})*Tp{nox+1}'*Up;    

%---------- extra stuff -----------------
% should work but not fully tested (MB 2007-02-19)
sstot_Y = sum(sum(Y.*Y));
F=Y-Up*Cp';
R2Y = 1 - sum(sum( F.*F))/sstot_Y;
%---------


EEprime=(K{nox+1,nox+1}-Tp{nox+1}*Tp{nox+1}');
sstot_K = (sum(diag(K{1,1})));
R2X = [];
R2XO = [];
R2XC = [];
R2Yhat = []; % R2Yhat 22 Jan 2010 / MR

for (i = 1:(nox+1))
    rss = sum(diag( K{i,i}-Tp{i}*Tp{i}' ));    
    R2X = [R2X, 1 - rss/sstot_K ];
    
    rssc = sum(diag( K{1,1}-Tp{i}*Tp{i}' ));    
    R2XC = [R2XC, 1 - rssc/sstot_K ];
    
    rsso = sum(diag( K{i,i} ));    
    R2XO = [R2XO, 1 - rsso/sstot_K ];
	
	% R2Yhat 22 Jan 2010 / MR - not fully tested
	Yhat = Tp{i} * Bt{i} * Cp';
	R2Yhat = [R2Yhat, 1 - sum(sum((Yhat-Y).^2))/sstot_Y];
end


%----------------------------------------



model.Cp=Cp;
model.Sp=Sp;
model.Up=Up;
model.Tp=Tp;
model.T=Tp{nox+1};
model.co=co;
model.so=so;
model.to=to;
model.toNorm=toNorm;
model.Bt=Bt;
model.A=A;
model.nox=nox;
model.K=K;

%convert struct to matrix
model.To=repmat(0, size(model.T,1),nox);
for (i = 1:length(to))
    model.To(:,i)=to{i};
end

%extra stuff
model.EEprime=EEprime;
model.sstot_K=sstot_K;
model.R2X=R2X;
model.R2XO=R2XO;
model.R2XC=R2XC;
model.sstot_Y=sstot_Y;
model.R2Y=R2Y;
model.R2Yhat=R2Yhat; % R2Yhat 22 Jan 2010 / MR


%Pre-processing
model.preProc.K=preProcK;
model.preProc.Y=preProcY;
model.preProc.paramsY=scaleParams;

model.class='kopls';


end
