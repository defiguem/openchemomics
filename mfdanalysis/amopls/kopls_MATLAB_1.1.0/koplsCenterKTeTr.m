function [KteTr]=koplsCenterKTeTr(KteTr,KtrTr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Centering function for the hybrid test/training kernel, which
% is constructed from the test matrix Xte and the training matrix
% Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
% (un-centered) training kernel to estimate mean values
% (see 'koplsKernel()' for details on constructing a kernel matrix).
%
% ** INPUT
% KteTr = Hybrid test/training kernel matrix;
%   KteTr = <phi(Xte), phi(Xtr)>.
% KtrTr = Training kernel matrix; Ktrain = <phi(Xtr), phi(Xtr)>.
%
% ** OUTPUT
% KteTr = The centered kernel matrix.
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


[ntr,ktr]=size(KtrTr);
[nte,kte]=size(KteTr);

Itrain=eye(ntr);
I_nTrain=ones(ntr,1);
nTrain=ntr;

I=eye(nte);
I_n=ones(nte,1);
n=nte;

KteTr = (KteTr-(1/nTrain)*I_n*I_nTrain' * KtrTr) * (Itrain-(1/nTrain).*I_nTrain*I_nTrain');
end
