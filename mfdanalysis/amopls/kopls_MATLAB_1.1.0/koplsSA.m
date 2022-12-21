function [Optset,settings]=koplsSA(Xtr,Ytr,A,oax,modelType,optargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform simulated annealing for optimisation of the kernel
% parameter. The error is mimimised, with the error as defined in
% koplsModelInternal, (1-Q2) mimimisation for regression and minimisation
% of (1-sensitivity) for discriminant analysis target function, or (1-area 
% under ROC curve) for 'daAUC' discriminant analysis. The correct 
% data set should be loaded with Xtr as X training set and Ytr for the
% Y training set, these data are used to perform a cross-validation
% on to determine the optimal kernel parameter setting. 
%
% Use:
% [Optset,settings]=koplsSA(Xtr,Ytr,A,oax,modelType,{optional, input})
%
% ** INPUT: 
% Required input parameters:
% Xtr = Training X data.
% Ytr = Training Y data.
% A = The number of Y-predictive components (integer).
% oax = The number of Y-orthogonal components (integer).
% modelType = 're' for regression, 'da' and 'daAUC' for discriminant 
%   analysis, if 'da', the mean sensitivity is optimised, for 'daAUC' the 
%   area under the receiver operating characteristic curve is optimised 
%   (only for two-class problems). 
%
% Possible input parameters for Simulated Annealing, pair-wise input in 
% optargin:
% kernelType = 'g' for gaussian (default).
% preProcK = Pre-processing settings of the X data.
%   Either 'mc' for mean-centering (default) or 'no' for no pre-processing.
% preProcY = Pre-processing settings for Y data. Either 'mc' for
%   mean-centering (default), 'uv' for mc + scaling to unit-variance,
%   'pa' for mc + Pareto-scaling or 'no' for no scaling.
% StX = Starting position for algorithm search - can be initiated
%   automatically or from a number of preset positions. 
% T_0 = Starting temperature; default = 0.1.
% epsilon = Termination criterion; default 0.01.
% Neps = Number of subsequent optimal points evaluated for convergence
%   criterion; default 2.
% v_0 = Starting step vector size, will be adjusted throughout optimisation 
%   but should preferable be able to cover a large range of the possible 
%   optimal kernel parameter values; default is StX.
% Ns = Number of points before reducing vector length; default 5.
% Nt = Number of vector reductions before temperature update; default 5.
% rT = Exponential cooling of temperature; default 0.1.
% nrreps = Number of runs; default 1.
% verbose = 1 for plotting and displaying temperature updates and results, 
%   0 for not (default).
%
% Possible input PLS model parameters, pair-wise input in optargin:
% nrcvinner = Number of cross-validation rounds (integer); default 10.
% cvType = Type of cross-validation. Either 'nfold' for n-fold
%   cross-validation, 'mccv' for Monte Carlo CV or 'mccvb' for Monte Carlo
%   class-balanced CV; default 'mccv' for 're' and 'mccvb' for 'da' and 
%   'daAUC' modelType.
% cvFrac = Fraction of observations in the trainingset during
%   crossvalidation. Only applicable to 'mccv', 'mccvb' crossvalidation;
%   default 0.75.
%
% ** OUTPUT:
% Optset = Optimal settings.
% settings = Settings and outcomes for simulated annealing.
%
% Reference: Corona et al; 'minimizing multimodal functions of continuous
% variables with the "simulated annealing" algorithm', ACM transactions on
% Mathematical Software, vol. 13, no.3, september 1987, pages 262-280.
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

%%%%% basic checks and default settings

if isempty(find(strcmp(optargin,'kernelType'))) 
    kernelType='g';
else 
    kernelType=optargin{(find(strcmp(optargin,'kernelType'))+1)};
end
if ~((strcmp(kernelType,'g')))
    disp('this is an incorrect setting for kernelType')
end

if isempty(find(strcmp(optargin,'preProcK'))) 
    preProcK='mc';
else 
    preProcK=optargin{(find(strcmp(optargin,'preProcK'))+1)};
end

if ~((strcmp(preProcK,'mc')) | (strcmp(preProcK,'no')))
    disp('this is an incorrect setting for X scaling')
end

if isempty(find(strcmp(optargin,'preProcY')))
    preProcY='mc';
else 
    preProcY=optargin{(find(strcmp(optargin,'preProcY'))+1)};
end
if ~((strcmp(preProcY,'uv')) | (strcmp(preProcY,'mc')) | (strcmp(preProcY,'no')) | (strcmp(preProcY,'pa')))
    disp('this is an incorrect setting for Y scaling')
end

if isempty(find(strcmp(optargin,'nrreps')))
    nrreps=1;
else 
    nrreps=optargin{(find(strcmp(optargin,'nrreps'))+1)};
end

if isempty(find(strcmp(optargin,'stX')))
    if strcmp(kernelType,'g')
        stX=rand(1); 
    else
        error('your kerneltype is incorrect, cannot determine best starting value')
    end
else
    stX=optargin{(find(strcmp(optargin,'stX'))+1)};
end
if ~(size(stX,2)==nrreps)
    error(['start values should have as many entries as number of repititions:  ' num2str(nrreps)])
elseif (min(stX)<=0)
    error(['start value should be bigger than zero'])
end

if isempty(find(strcmp(optargin,'T_0')))    
    T_0=0.1;
else
    T_0=optargin{(find(strcmp(optargin,'T_0'))+1)};
end

if isempty(find(strcmp(optargin,'epsilon')))    
    epsilon=0.01;
else
    epsilon=optargin{(find(strcmp(optargin,'epsilon'))+1)};
end

if isempty(find(strcmp(optargin,'Neps')))    
    Neps=2;
else
    Neps=optargin{(find(strcmp(optargin,'Neps'))+1)};
end

if isempty(find(strcmp(optargin,'v_0')))    
    v_0=stX;
else
    v_0=optargin{(find(strcmp(optargin,'v_0'))+1)};
end

if isempty(find(strcmp(optargin,'Ns')))    
    Ns=5;
else
    Ns=optargin{(find(strcmp(optargin,'Ns'))+1)};
end

if isempty(find(strcmp(optargin,'Nt')))    
    Nt=5;
else
    Nt=optargin{(find(strcmp(optargin,'Nt'))+1)};
end

if isempty(find(strcmp(optargin,'rT')))    
    rT=0.1;
else
    rT=optargin{(find(strcmp(optargin,'rT'))+1)};
end
if (rT>1 | rT<0)
    error('the rate of decay in temperature rT should be between 0 and 1')
end

if isempty(find(strcmp(optargin,'nrcvinner')))    
    nrcvinner=10;
else
    nrcvinner=optargin{(find(strcmp(optargin,'nrcvinner'))+1)};
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

if isempty(find(strcmp(optargin,'verbose')))
    verbose=0;
else
    verbose=optargin{(find(strcmp(optargin,'verbose'))+1)};
end

if (isempty(modelType) | ((~strcmp(modelType,'daAUC'))&(~strcmp(modelType,'da')) & (~strcmp(modelType,'re'))))
    error('this is not a valid modelType')
end

if isempty(Xtr)
    error('this is not a valid X data matrix')
end

if isempty(Ytr)
    error('this is not a valid Y data matrix')
end

if isempty(A)
    error('this is not a correct setting for the number of Y-predictive components, A, (integer).')
end

if isempty(oax)
    error('this is not a correct setting for the number of Y-orthogonal components, oax, (integer).')
end


%%%%% Start of the iterations

for reps=1:nrreps
%%%%% Step 0: Initialisation
timea=cputime;
f_0=koplsModelInternal(Xtr,Ytr,A,oax,nrcvinner,cvType,preProcK,preProcY,cvFrac,modelType,kernelType,stX(nrreps));

f_old=f_0;
X_old=stX(nrreps);
v_new=v_0;
n_1=0;
T_new=T_0;
c_1=2; 
maxits=1000; %maximum number of iterations to prevent long calculation times.
Xlist=[X_old];%list of all evaluated Xs
flist=[f_old];%list of all evaluated fs
Tlist=[T_new];%list of all evaluated temperatures
X_opt=X_old;
f_opt=f_old;
Xoptlist=[X_opt];%list of increasingly better Xs
foptlist=[f_opt];%list of increasingly better fs
Xshortlist=[X_opt];%list of Xs after each temperature change
fshortlist=[f_opt];%list of fs after each temperature change
vlist=[v_new]; %list of step vectors
Tshortlist=[T_new]; %list of Temperatures.


%%%%% Step 1: Calculate new X
teli=1; 
telk=1; 

while (not((length(fshortlist)>Neps) && (min((fshortlist(end-Neps:end-1))-fshortlist(end)<=epsilon) & (fshortlist(end)-f_opt)<=epsilon))) 
telm=1;
while telm<=Nt % time to reduce the temperature otherwise
    telj=1; 
    while telj<=Ns    %total number of points before reducing vector length
        X_new=X_old+(2*rand(1)-1)*v_new;    %starting from point X_old, generate a random point Xnew
        while X_new<=0 ; %if new point is negative, try again
            X_new=X_old+(2*rand(1)-1)*v_new; 
        end

        %%%%% Step 2: Evaluate value in new X, accept or reject this as new starting point
        f_new=koplsModelInternal(Xtr,Ytr,A,oax,nrcvinner,cvType,preProcK,preProcY,cvFrac,modelType,kernelType,X_new); %Evaluate function for new X
        if f_new <=f_old % If improved, update data tables
            X_old=X_new;
            f_old=f_new;
            Xlist=[Xlist,X_old];
            flist=[flist,f_old];
            Tlist=[Tlist,T_new];
            teli=teli+1;
            n_1=n_1+1; %used to determine how to update step vector; number of accepted moves  
            if f_old<f_opt 
                X_opt=X_old; 
                f_opt=f_old;
                Xoptlist=[Xoptlist,X_opt]; 
                foptlist=[foptlist,f_opt];
            end

        else %If not improved, accept or reject slightly worse point with probability from Metropolis move
            compnr=rand(1);
            if compnr<exp((f_old-f_new)/T_new)
                X_old=X_new; 
                f_old=f_new;
                Xlist=[Xlist,X_old];
                flist=[flist,f_old];
                Tlist=[Tlist,T_new];
                teli=teli+1;
                n_1=n_1+1;
            end
        end

        telj=telj+1; 
    end

    %%%%% Step 3: Change the step vector component, 
    if n_1>0.6*Ns 
        v_new=v_new*(1+c_1*(((n_1/Ns)-0.6)/0.4));
    elseif n_1 < 0.4*Ns
        v_new=v_new/(1+c_1*((0.4-(n_1/Ns))/0.4));
    else
        v_new=v_new;
    end

    vlist=[vlist,v_new]; 
    n_1=0;
    telm=telm+1;
end

%%%%% Step 4: Update temperature
T_new=rT*T_new; 
Tshortlist=[Tshortlist,T_new];
if verbose==1
    disp(['Temperature update: ' num2str(T_new)])
end
Xshortlist=[Xshortlist,X_old]; 
fshortlist=[fshortlist,f_old];
if verbose==1
    figure(345);subplot(2,1,1),plot(Tshortlist,'bx-');title('Temperature','Fontsize',16);subplot(2,1,2);plot(fshortlist,'rx-');title('Optimization','Fontsize',16);xlabel('number of temperature updates')
end
X_old=X_opt;
f_old=f_opt; 
telk=telk+1;
teli=teli+1;

%%%%% Step 5: stop after convergence or exceeding maximum number of
%%%%% iterations
if teli>maxits
    disp('stopped before convergence')
    break
end
end


Optset(reps)=X_opt(end); %optimal settings
settings.time{reps}=cputime-timea; % total time per optimisation
settings.Q2{reps}=1-f_opt(end); % Q2 or sensitivity if DA per optimisation
settings.preProcK=preProcK; % scaling of X matrix
settings.preProcY=preProcY; % scaling of Y matrix
settings.stX(reps,:)=stX(reps); % starting value of kernel parameter
settings.T_0=T_0; % starting temperature
settings.epsilon=epsilon; %Termination criterion
settings.Neps=Neps; %Number of subsequent optimised points evaluated for convergence
settings.v_0(reps,:)=v_0; %Starting step vector size 
settings.Ns=Ns; %Number of points before reducing vector length
settings.Nt=Nt; %Number of vector reductions before temperature update
settings.rT=rT; %rate of exponential cooling of temperature
settings.nrreps=nrreps; %Number of runs
settings.nrcvinner=nrcvinner; % number of cross-validations
settings.cvType=cvType; %type of cross-validation
settings.cvFrac=cvFrac; %Fraction of observations in the trainingset during cross-validation
settings.modelType=modelType; %PLS model type: regression or discriminant analysis
settings.Xtr=Xtr; %training X data
settings.Ytr=Ytr; %training Y data
settings.A=A; % number of predictive components
settings.oax=oax; % number of Y-orthogonal components
settings.n_1=n_1; %number of optimisation directions
settings.c_1=c_1; % varying criterion for step size update
settings.kernelType=kernelType; %kernel type
% settings.nrers=nrers; % number of convergence warnings
settings.maxits=maxits; %maximum number of iterations
settings.Xlist{reps}=Xlist; %list of all evaluated Xs
settings.flist{reps}=flist; %list of all evaluated fs
settings.Tlist{reps}=Tlist;%JMF added 25 may 2010
settings.Xoptlist{reps}=Xoptlist; %list of increasingly better Xs
settings.foptlist{reps}=foptlist;% list of increasingly better fs
settings.Xshortlist{reps}=Xshortlist; %list of Xs after each temperature change
settings.fshortlist{reps}=fshortlist; %list of fs after each temperature change
settings.Tshortlist{reps}=Tshortlist; %list of temperature update
settings.vlist{reps}=vlist; %list of step vectors
end
end



