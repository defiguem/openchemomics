function [sensvec, specvec, classvec, tot_sens,meanSens,meanSpec]=koplsSensSpec(v, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates sensitivity and specificity in a class-wise fashion.
%
% ** INPUT
% v = row vector of true class assignments (template).
% m = matrix (or row vector) of class assignments to be compared.
%
% ** OUTPUT
% sensvec = sensitivity for each class.
% specvec = specificity for each class.
% classvec = the class identifier corresponding to each column in
%   sensvec and specvec.
% tot_sens = total sensitivity.
% meanSens = mean sensitivity over all classes.
% meanSpec = mean specificity over all classes.
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


%Make sure data seems OK
[nrow, ncol] = size(m);
if (nrow ~= length(v))
	disp('Template vector differs in length from the matrix to be compared. Bailing out.')
	return
end


classes = sort(unique(v));


%append multiple columns to the end of v and m
if (ncol > 1)
	
	%disp('More than 1 column in m... ')
	
	mtemp = [];
	vtemp = [];
	for (j = 1:nrow)
		for (i = 1:ncol)
		
			%the first column is always added
			if (i == 1)
				vtemp = [vtemp; v(j)];
				%if we have NaN values, we claim it is a class
				%that does not exist --> incorrect pred.
				if (isnan(m(j, i)))
					%disp('Adding value..')
					mtemp = [mtemp; max(classes)+1];
					
				else
					mtemp = [mtemp; m(j,i)];
				end
						
			else
				if (~isnan(m(j, i)))
					vtemp = [vtemp; v(j)];
					mtemp = [mtemp; m(j, i)];
				end
			end
			
		end
	end
	
	v=vtemp;
	m=mtemp;
end

row_indices = 1:1:length(v);

tot_sens = sum( v == m)/length(v);


ind=find(~isnan(m));
%%'number of observations that have not been assigned a class:';
%%sum(isnan(m));
m=m(ind);
v=v(ind);
for (i = 1:length(classes))
	
	class = classes(i);
	rows_inclass = row_indices( v == class);
	rows_notinclass = row_indices( v ~= class);
	
	tp_rate = sum( m(rows_inclass) == class) / length(rows_inclass);
	tn_rate = sum( m(rows_notinclass) ~= class) / length(rows_notinclass);


	sensvec(i) = tp_rate;
	specvec(i) = tn_rate;
	classvec(i) = class;
	
end


meanSens=mean(sensvec);
meanSpec=mean(specvec);
