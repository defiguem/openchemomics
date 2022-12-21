function [dummy, labels_sorted]=koplsDummy(class, numClasses);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Converts integer vector to binary matrix (dummy matrix).
%
% ** INPUT
% class = vector with class belongings (integer).
% numClasses = pre-defines the number of classes in the output
%	(if undefined, the number of unique entries in 'class' will be
%   used).
%
% ** OUTPUT
% dummy = A matrix with rows corresponding to observations and
%   columns to classes. Each element in matrix is either one
%   (observation belongs to class) or zero (observation does not
%   belong to class).
% labels_sorted = The class labels that are found in class in
%	sorted order.
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

labels=unique(class);%the set of classlabels in class
labels_sorted=sort(labels); %sort labels in ascending order

len_class=length(class);%number of samples

%number of classes
if (nargin < 2)
	len_labels=length(labels);
else
	len_labels=numClasses; %%force number of classes
end

dummy=zeros(len_class,len_labels); %dummy matrix initialized as a zero matrix

for i=1:len_labels %for each class label
   ind=find(class==labels_sorted(i)); %find the rows (samples) that belongs to the current class, labels_sorted(i
   dummy(ind,i)=1; %write ones in the positions where we have the current class....
end


