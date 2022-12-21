function []=koplsPlotModelDiagnostics(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots diagnostic parameters from K-OPLS cross-validation.
% This includes:
% - R2X = Cumulative explained variation for all model components.
% - R2XO = Cumulative explained variation for all Y-orthogonal model
%   components.
% - R2XC = Explained variation for predictive model components after
%       addition of Y-orthogonal model components.
%
% For further information regarding the definitation and calculation
% of these quantities, see e.g appendix of:
%  * Trygg J and Wold S. J Chemometrics 2003; 17:53-64.
%
% ** INPUT
% model = a model constructed using 'koplsModel()'.
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

if (~strcmp(model.class, 'kopls'))
    error('Model must be of type "kopls". Aborting.')
end

xlabs=cell(1);
for (i = 1:model.A)
    xlabs{1}=strcat('Tp,',num2str(i));
end

xlabsOrtho=cell(1);
if (model.nox > 0)
    for (i = 1:size(model.To,1))
        xlabsOrtho{i}=strcat('To,',num2str(i));
        xlabs{i+model.A}=xlabsOrtho{i};
    end
end

subplots.x=2;
subplots.y=2;

%%Plot
subplot(subplots.x,subplots.y,1);
bar(model.R2X, 'g')
ylabel('R2X (cumulative)');
set(gca,'XTickLabel',xlabs)
%, names.arg=labels, ylab="R2X (cumulative)", col="Green")

subplot(subplots.x,subplots.y,2);
bar(model.R2XO(2:length(model.R2XO)), 'b')
ylabel('R2Xortho (cumulative)');
set(gca,'XTickLabel',xlabsOrtho)
%, names.arg=labels[2:length(model.R2XO)], ylab="R2Xortho (cumulative)", col="Blue")

subplot(subplots.x,subplots.y,3);
bar(model.R2XC, 'r');
ylabel('R2Xcorr');
set(gca,'XTickLabel',xlabs)
%, names.arg=labels, ylab="R2Xcorr", col="Red")

% Evaluate if Q2 can be plotted, if not display an 'error plot'
subplot(subplots.x,subplots.y,4);
try
    bar(model.Q2', 'y');
    ylabel('Q2Y');
    set(gca,'XTickLabel',xlabs)
    %, names.arg=labels, ylab="Q2", col="Yellow")
catch
    text(0.1, 0.5, 'Q2Y undefined (cross-validation not performed)', 'FontSize',8);
end

% create main header
ax=axes('Position',[.08 .08 .855 .855],'Visible','off');
set(get(ax,'Title'),'Visible','on');
title('K-OPLS model diagnostics');
  
%subplot(subplots.x, subplots.y, 1:4);


