function [scaleS]=koplsScale(X,centerType,scaleType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for mean-centering and scaling of a matrix.
% 
% ** INPUT
% X = X matrix (to be mean-centered/scaled).
% centerType = 'mc' for mean-centering, 'no' for no centering.
% scaleType = 'uv' for unit variance scaling, 'pa' for Pareto
%   scaling, 'no' for no scaling.
%
% ** OUTPUT:
% scaleS = An object with the following properties:
%   centerType = 'mc' or 'no'.
%   scaleType = 'uv', 'pa' or 'no'.
%   meanV = vector with mean values for all columns in X.
%   stdV = vector with standard deviations for all columns in X.
%   X = Original input matrix X, scaled according to 'centerType'
%       and 'scaleType'.
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


    scaleS.centerType=centerType;
	scaleS.scaleType=scaleType;
	scaleS.meanV=mean(X);
	scaleS.stdV=std(X);
	[m,n]=size(X);	
	if(centerType=='mc')
		X=X-ones(m,1)*scaleS.meanV;
	end
	if(scaleType=='uv')
		X=[X./repmat(scaleS.stdV,m,1)];
    end
    
	if(scaleType=='pa')
		X=[X./repmat(sqrt(scaleS.stdV),m,1)];
    end	
	scaleS.X = X;
	return;	
end





