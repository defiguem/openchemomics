function [res]=koplsPlotSensSpec(modelFull)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots sensitivity and specificity results from cross-validation
% in a bar plot. The produced bars are shown separately for each
% class including overall sensitivity and specificity results.
%
% ** INPUT
% modelFull = the K-OPLS cross-validation results from 'koplsCV()'
%
% ** OUTPUT
% res = The resulting sensitivity and specificity measures.
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


if (~strcmp(modelFull.class, 'koplscv'))
    error('Unknown model type (must be of type "koplscv"). Aborting.');
end

labels=[];
try
	num_entries=length(modelFull.da.meanSensAllOsc);
	plot_mat=repmat(0, 2, num_entries );
	
	for (i = 1:num_entries)
		plot_mat(1, i)=modelFull.da.meanSensAllOsc{i};
		plot_mat(2, i)=modelFull.da.meanSpecAllOsc{i};
		labels=[labels; strcat('To,', num2str(i-1) )];
	end
	
    %[sensvec, specvec, classvec, tot_sens]=koplsSensSpec(modelFull.da.trueClass, modelFull.da.predClass);
	

catch
    error('The cross-validation results are not for discriminant analysis. Aborting');
end
    
h=bar(plot_mat');
h=gca;
set(h, 'XTickLabel', labels)
ylabel('Sens. and spec. (%)');

sensvec=plot_mat(1,:);
specvec=plot_mat(2,:);
res.sens=sensvec;
res.spec=specvec;
res.tot_sens=modelFull.da.tot_sens;
res.meanSens=mean(sensvec);
res.meanSpec=mean(specvec);

% create main header
ax=axes('Position',[.08 .08 .855 .855],'Visible','off');
set(get(ax,'Title'),'Visible','on');
title('Sensitivity and specificity over cross-validation');

