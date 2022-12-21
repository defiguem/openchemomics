function []=koplsPlotCVDiagnostics(modelFull)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots diagnostic parameters from K-OPLS cross-validation.
% This includes:
% - R2X = Cumulative explained variation for all model components.
% - R2Xortho = Cumulative explained variation for all Y-orthogonal
% 	model components.
% - R2Xcorr = Explained variation for predictive model components
% 	after addition of Y-orthogonal model components.
% - Q2Y = Total Q-square result (1 - pred. residual / original var.)
%   for all Y-orthogonal components.
%
% For further information regarding the definitation and calculation
% of these quantities, see e.g the appendix of:
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

if (~strcmp(modelFull.class, 'koplscv'))
    error('Unknown model type (must be of type "koplscv"). Aborting.')
end

model=modelFull.koplsModel;
model.Q2=modelFull.cv.Q2Yhat;

koplsPlotModelDiagnostics(model)



