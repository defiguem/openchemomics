function [scaleS]=koplsRescale(scaleS,varargin)	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scales a matrix based on pre-defined parameters from a scaling
% object.
%
% ** INPUT
% scaleS = An object containing scaling parameters
%   (see 'koplsScale()').
% varargin = If defined, this matrix will be scaled and returned.
%	Otherwise the original data set in the scaleS object will be
%   scaled and returned.
%
% ** OUTPUT
% scaleS = An object containing the following entries:
%   centerType = 'mc' (mean-centering) or 'no' (no centering).
%   scaleType = 'uv' (unit variance), 'pa' (pareto) or 'no'
%       (no scaling).
%   meanV = vector with mean values for all columns in X.
%   stdV = vector with standard deviations for all columns in X.
%   X = Scaled version of 'varargin', if defined, otherwise,
%       scaled version of scaleS.X from input. Scaling is done
%       according to 'centerType' and 'scaleType'.
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


	if(~isempty(varargin))
        X=varargin{1};
    else
        X=scaleS.X;
    end
    [m,n]=size(X);
	
	
	if(scaleS.scaleType=='uv')
		X=X.*repmat(scaleS.stdV,m,1);
	end
	if(scaleS.scaleType=='pa')
		X=X.*repmat(sqrt(scaleS.stdV),m,1);
	end
	
	if(scaleS.centerType=='mc')
		X=X+ones(m,1)*scaleS.meanV;
	end
	
	scaleS.centerType='no';
	scaleS.scaleType='no';
	scaleS.X=X;
	return;	
end



