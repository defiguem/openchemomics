function [scaleSA]=koplsScaleApply(X,scaleS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applies scaling from external scaling objects (see
% 'koplsScale()') on a matrix X. Returns the scaled matrix,
% including scaling settings.
% 
% ** INPUT
% X = X matrix (to be scaled).
% scaleS = An object containing scaling parameters
%   (see 'koplsScale()').
% 
% ** OUTPUT
% scaleSA = An object with the following properties:
%   centerType = 'mc' or 'no'.
%   scaleType = 'uv', 'pa' or 'no'.
%   meanV = vector with mean values for all columns in X.
%   stdV = vector with standard deviations for all columns in X.
%   X = Scaled version of 'X', scaled according to 'centerType'
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


	scaleSA.centerType=scaleS.centerType;
	scaleSA.scaleType=scaleS.scaleType;
	scaleSA.meanV=scaleS.meanV;
	scaleSA.stdV=scaleS.stdV;
	[m,n]=size(X);
	
	if(scaleS.centerType=='mc')
		X=X-ones(m,1)*scaleS.meanV;
	end
	if(scaleS.scaleType=='uv')
		X=X./repmat(scaleS.stdV,m,1);
	end
	if(scaleS.scaleType=='pa')
		X=[X./repmat(sqrt(scaleS.stdV),m,1)];
	end
	
	scaleSA.X=X;
	return;	
end
