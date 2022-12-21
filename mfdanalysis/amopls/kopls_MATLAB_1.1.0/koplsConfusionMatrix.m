function [A]=koplsConfusionMatrix(true,pred)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates a confusion matrix from classification results.
%
% ** INPUT
% true = true class belonging.
% pred = predicted class assignment.
%
% ** OUTPUT
% A = Confusion matrix.
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


uniqueClass=unique(true);

A=zeros(length(uniqueClass),length(uniqueClass));

for(i=1:length(uniqueClass))%for each class
	indTrue=find(true==uniqueClass(i));
		for(j = 1:length(indTrue))
			%pred(indTrue(j))
			%[find(uniqueClass==pred(indTrue(j)))\
			A(i,find(uniqueClass==pred(indTrue(j))))=A(i,find(uniqueClass==pred(indTrue(j))))+1;
        end
        A(i,:)=A(i,:)./length(indTrue);
end

end
