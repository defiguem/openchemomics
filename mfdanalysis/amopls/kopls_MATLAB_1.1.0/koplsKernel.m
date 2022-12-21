function [K]=koplsKernel(X1,X2,Ktype,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constructs a kernel matrix K = <phi(X1), phi(X2)>.
% The kernel function k() determines how the data is transformed and
% is passed as the separate parameter 'Ktype' to the function.
% Currently 'Ktype' can be either 'g' (Gaussian) or 'p' (polynomial).
%
% ** INPUT
% X1 = the first X matrix (non-centered). This is the left side in
%   the expression K = <phi(X1), phi(X2)>.
% X2 = the second X matrix (non-centered). This is the right side
%   in the expression K = <phi(X1), phi(X2)>.
%  If X2 = [] (empty set), then only X1 will be used for
%  the calculations. This way, only (n^2 - n)/2 instead of n^2
%  calculations have to be performed, which is typically much
%  faster. Only applicable for pure training or testing kernels.
% Ktype = the type of kernel used. Supported entries are:
%   - 'g': Gaussian kernel.
%   - 'p': Polynomial kernel.
% params = A vector with parameter for the kernel function.
%   (Currently, all supported kernel functions use a scalar value
%    so the vector property of the parameters is for future
%    compability).
%
% ** OUTPUT
% K = The kernel matrix, transformed by the kernel function
%   specified by 'Ktype'.
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



    [l1, r1] = size(X1);
    [l2, r2] = size(X2);

	fast_calc = (l2 == 0); %% we only have one matrix (i.e. X1 == X2)
	
	if (fast_calc)
		K = repmat(0, l1, l1);
	else
		K = repmat(0, l1, l2);
	end

    if(Ktype=='g') %gaussian
        sigma=params(1); %small value = overfit, larger = more general

        if (fast_calc) 
            
            % calculate upper triangle and duplicate (due to symmetry)
            for(i = 1:l1)
                for( j = i:l1)
					K(i,j)=exp(- (norm( X1(i,:) - X1(j,:) )^2 ) / (2*sigma^2) );
                    K(j,i)=K(i,j);
                    %%wikipedia ntation - same -same, just different parameter
                    %%(sigma...)
                    %K(i,j)=exp(- (norm( X1(i,:) - X2(j,:) )^2 ) / sigma ); %rosipal match
                end

            end

        else
             % loop over entire kernel matrix
             for(i = 1:l1)
                for( j = 1:l2)
                    K(i,j)=exp(- (norm( X1(i,:) - X2(j,:) )^2 ) / (2*sigma^2) );
                    %%wikipedia ntation - same -same, just different parameter
                    %%(sigma...)
                    %K(i,j)=exp(- (norm( X1(i,:) - X2(j,:) )^2 ) / sigma ); %rosipal match
                end

            end            
        end
        
    end


    if(Ktype=='p') %polynomial, order=param
        %sigma=0.0005; %small value = overfit, larger = more general
        porder=params(1);
        
        
        if (fast_calc) %% we only have one matrix (i.e. X1 == X2)
            
             % calculate upper triangle and duplicate (due to symmetry)
            for(i = 1:l1)
                for( j = i:l1)		
                    K(i,j)=(X1(i,:) * X1(j,:)' +1 ).^porder; %polynomial/homogenous
                    K(j,i)=K(i,j);
                    %K(i,j)=(X1(i,:) * X2(j,:)' -1)^porder; %polynomial non-homogenous
                end
            end
        else

            % loop over entire kernel matrix
            for(i = 1:l1)
                for( j = 1:l2)		
                    K(i,j)=(X1(i,:) * X2(j,:)' +1 ).^porder; %polynomial/homogenous
                    %K(i,j)=(X1(i,:) * X2(j,:)' -1)^porder; %polynomial non-homogenous
                end
            end
        end
    end

end
