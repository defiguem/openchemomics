function []=koplsPlotScores(model,x,xsub,y,ysub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Produces score plots from K-OPLS models. If model components are
% unspecified, all possible combinations are displayed as a scatter
% plot matrix (default). Otherwise, two selected components will be
% shown using a traditional 2D scatter plot.
%
% NB: The diagonal of the scatter plot matrix requires functionality
% from the 'MATLAB statistics toolbox' to produce variable densities.
% In the absence of the toolbox, the diagonal will consist of blank
% plots.
%
% ** INPUT
% model = the K-OPLS model generated using 'koplsModel()'.
% x = the vector index for the x axis.
% xsub = the vector identifier {'p', 'o'} for the x axis.
% y = the vector index for the y axis.
% ysub = the vector identifier {'p', 'o'} for the y axis
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

plotAll=(nargin < 2);

% The default visualization is a scatter plot matrix
% of small multiples.
% Diagonal depicts the density of a variable.
if (plotAll)

    Tall=model.T;
    xlabs=[];
    for (i = 1:size(Tall,2))
        xlabs=[xlabs; strcat('Tp,', num2str(i))];
    end
            
    %labels=paste("tp", 1:model.A, sep=",")

    if (model.nox > 0)
        Tall=[Tall, model.To];
        for (i = 1:model.nox)
            xlabs=[xlabs; strcat('To,', num2str(i))];
        end
        %labels=c(labels, paste("to", 1:model.nox, sep=","))
    end
    
    subplots.x=size(Tall,2);
    subplots.y=size(Tall,2);

    ind=0;
    for (i = 1:size(Tall,2))
        for (j = 1:size(Tall,2))

            %keyboard;

            ind=ind+1;
            subplot(subplots.x,subplots.y,ind);

            if (i == j)
                try
					[f,xi]=ksdensity(Tall(:,i));
					plot(xi,f);
					xlabel(xlabs(i,:))
					ylabel('Density')
					%plot(density(Tall[, i]), main="", ylab="Density", xlab=labels[i])
				catch
					%disp('ksdensity() failed')
					text(0.1, 0.5, 'Stat. toolbox required to see density', 'FontSize',8);
				end
            else
                plot(Tall(:,i), Tall(:,j), 'xb')
                xlabel(xlabs(i,:))
                ylabel(xlabs(j,:))
            end
        end
    end
    
    % create main header
    ax=axes('Position',[.08 .08 .855 .855],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    title('Scatter plot matrix of K-OPLS model scores');

else

    if ( (strcmp(xsub,'p') && x > model.A) || (strcmp(xsub,'o') && x > model.nox) )
        error('X variable outside range of model. Aborting.')
    end

    if ( (strcmp(ysub,'p') && y > model.A) || (strcmp(ysub,'o') && y > model.nox) )
        error('Y variable outside range of model. Aborting.')
    end

    if (strcmp(xsub,'p'))
        xvec=model.T(:,x);
    elseif (strcmp(xsub,'o'))
        xvec=model.To(:,x);
    else
        error('Unknown model component specification for x: should be {o, p}. Aborting.')
    end

    if (strcmp(ysub,'p'))
        yvec=model.T(:,y);
    elseif (strcmp(ysub,'o'))
        yvec=model.To(:,y);
    else
        error('Unknown model component specification for y: should be {o, p}. Aborting.')
    end

    if (x==y && xsub==ysub)
        try
			[f,xi]=ksdensity(xvec);
			plot(xi,f);
			xlabel(strcat('T',xsub,',',num2str(x)))
			ylabel('Density');     
		catch
			%disp('ksdensity() failed')
			text(0.1, 0.5, 'Stat. toolbox required to see density', 'FontSize',8);
		end
    else
        plot(xvec, yvec, 'xb')
        xlabel(strcat('T',xsub,',',num2str(x)))
        ylabel(strcat('T',ysub,',',num2str(y)))
    end

    % create main header
    ax=axes('Position',[.08 .08 .855 .855],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    title('Scatter plot of K-OPLS model scores');
    
end




