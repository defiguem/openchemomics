function [hh]=koplsPlotOptResults(model,optmethod,modelType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for plotting the results of the nested cross-validation for
% optimisation of the kerenel parameter. Plotted are Q2Yhat and the R
% values as a function of the number of orthogonal components for the final
% model after optimisation. In the right-hand subplot, the kernel parameter
% optimisation for the different CV rounds. Note that only temporary best 
% values are plotted, not the whole of the searched space for simulated
% annealing. For simulated annealing, a bar plot is created demonstrating 
% the time the individual CV rounds took.
%
% ** INPUT
% model = the returned cross-validated model (gridsearch or simulated
%   annealing)
% optmethod = 'SA' or 'GS', method used for optimisation
% modeltype = 're' for regression, 'da' for discriminant analysis optimised
%   to (mean sensitivity), 'daAUC' for optimised area under receiver
%   operating characteristic curve.
%
% ** OUTPUT
% [hh] = handle for produced figure
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
nroc=length(model.cv.Q2Yhat)-1;

hh=figure;
subplot(2,3,1)
hold on
plot(0:nroc,model.cv.Q2Yhat,'bx-')
plot(0:nroc,model.koplsModel.R2X,'r+-')
plot(0:nroc,model.koplsModel.R2XC,'m*-')
plot(0:nroc,model.koplsModel.R2XO,'go-')
legend('Q2Yhat','R2X','R2XC','R2XO','Location','Best')
xlabel('# orthogonal components')
title(['final model parameters'])
axis tight

subplot(2,3,[2,3,5,6]);
hold on
if strcmp(optmethod,'SA')
    count1=length(model.SAsettings);
    colmap=jet(count1);
    for teli=1:count1
        plot(model.SAsettings{1,teli}.Xoptlist{1,1},model.SAsettings{1,teli}.foptlist{1,1},'x-','Color',colmap(teli,:))
        legendstr(teli)={['CV' num2str(teli)]};
    end
    plot([model.KParamfinal model.KParamfinal],[0 1],'k:','Linewidth',3)
    legend([legendstr,['Full model:' num2str(model.KParamfinal,3)]],'Location','Best')
    for teli=1:count1
        plot(model.SAsettings{1,teli}.Xoptlist{1,1}(end),model.SAsettings{1,teli}.foptlist{1,1}(end),'.','Color',colmap(teli,:),'Markersize',30)
        timevec(teli)=model.SAsettings{1,teli}.time{1};    
    end
    xlabel('Kernel Parameter')
    if strcmp(modelType,'re')
        ylabel('1-Q^2_Y')
    elseif strcmp(modelType,'da')
        ylabel('1-(mean sensitivity)')
    elseif strcmp(modelType,'daAUC')
        ylabel('1-(AUROC)')
    end

    title('Optimisation results for individual cross-validation rounds for simulated annealing')
    axis tight
    
    subplot(2,3,4);
    bar(1:count1,timevec,'k')
    xlabel('CV round')
    ylabel('time (s)')
    title('time per CV')
    axis tight
end

if strcmp(optmethod,'GS')
    count1=length(model.GSsettings);
    colmap=jet(count1);
    for teli=1:count1
        plot(model.GSsettings{teli},model.GSresults{teli},'x-','Color',colmap(teli,:))
        legendstr(teli)={['CV' num2str(teli)]};
    end
    plot([model.KParamfinal model.KParamfinal],[0 1],'k:','Linewidth',3)
    legend([legendstr,['Full model:' num2str(model.KParamfinal,3)]],'Location','SouthWest')
    xlabel('Kernel Parameter')
    if strcmp(modelType,'re')
        ylabel('1-Q^2_Y')
    elseif strcmp(modelType,'da')
        ylabel('1-(mean sensitivity')
    elseif strcmp(modelType,'daAUC')
        ylabel('1-(AUROC)')
    end
    title('Optimisation results for individual cross-validation rounds for gridsearch')
    axis tight
end






