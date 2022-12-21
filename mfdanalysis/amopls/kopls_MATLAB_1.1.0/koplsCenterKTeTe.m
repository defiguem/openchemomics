function [KteTe]=koplsCenterKTeTe(KteTe,KteTr,KtrTr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Centering function for the test kernel, which is constructed
% from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
% Requires additional (un-centered) kernels KteTr and KteTr to
% estimate mean values (see 'koplsKernel()' for details on
% constructing a kernel matrix).
%
% ** INPUT
% KteTe = Test kernel matrix; KteTe = <phi(Xte), phi(Xte)>.
% KteTr = Test/training kernel matrix;
%   KteTr = <phi(Xte), phi(Xtr)>.
% KtrTr = Training kernel matrix; KtrTr = <phi(Xtr), phi(Xtr)>.
%
% ** OUTPUT
% KteTe = The centered test kernel matrix.
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


[nte,ntr]=size(KteTr);

Itrain=eye(ntr);
I_nTrain=ones(ntr,1);
nTrain=ntr;

I=eye(nte);
I_n=ones(nte,1);
n=nte;

D_te = (1/nTrain).*I_n*I_nTrain';
KteTe = KteTe - D_te*KteTr' - KteTr*D_te' + D_te*KtrTr*D_te';
end
