function [ROC]=koplsRoc(data,trueclass)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Calculates the approximate area under the receiver operating
% characteristic curve to assess the classification performance.
%
% ** INPUT
% data = matrix containing the predicted response matrix Y,
%   where columns denote classes and rows observations.
% trueclass = true class designation.
%
% ** OUTPUT
% ROC = the sensitivities and 1-specificities for different thresholds
%       as well as the area under curve
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
if (size(data,2))>2|(length(unique(trueclass))>2)
    error('ROC only possible for two-class problems')
end

classind=unique(trueclass);
ind1=find(trueclass==classind(1));
ind2=find(trueclass==classind(2));
tmp=data(:,1);
[a,b]=sort(tmp);
counter=1;
for trs=[a'];
    pred=2*ones(length(tmp),1);
    bigind=find(tmp>trs);
    pred(bigind)=1;  
    tp=length(find(pred(ind1)==1));        
    tn=length(find(pred(ind2)==2));
    fn=length(find(pred(ind1)==2));
    fp=length(find(pred(ind2)==1));
    sens(counter)=tp/(tp+fn);
    spec(counter)=1-(tn/(tn+fp));
    counter=counter+1;
end 
ROC.sens=sens;
ROC.spec=spec;
ROC.AUC=sum((spec(1:end-1)-spec(2:end)).*(sens(1:end-1)+sens(2:end))/2);

