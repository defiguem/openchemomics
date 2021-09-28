function [RV2] = rvmod(X,Y)
%
% The modified RV-coefficient can be used to characterize the relationship
% between pairs of high-dimensional matrices through a single number. This
% matrix correlation metric can be used in the same way as Pearson's
% correlation.
% 
% References :
% ============
% Smilde, A. K., Kiers, H. A. L., Bijlsma, S., Rubingh, C. M., & van Erk, 
% M. J. (2009). Matrix correlations for high-dimensional data : The 
% modified RV-coefficient. Bioinformatics, 25(3), 401‑405. 
% https://doi.org/10.1093/bioinformatics/btn634
% 
% Check section 2.3 of the above reference for details on calculation.
%
% Input arguments :
% =================
% X and Y : two high-dimensional matrices with the same number of 
%   observations in the rows and variables in the columns;
% 
% Output arguments :
% ==================
% RV2 : the modified RV-coefficient
% 
% Usage :
% =======
% [RV2] = rvmod(X,Y);
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
% 
% Modifications:
% ==============
%
% =========================================================================

XX = X*X';
YY = Y*Y';

XXtilde = XX - diag(diag(XX),0);
YYtilde = YY - diag(diag(YY),0);

vecXX = XXtilde(:);
vecYY = YYtilde(:);

num = vecXX' * vecYY;
denom = sqrt( (vecXX' * vecXX) * (vecYY' * vecYY) );
RV2 =  num / denom; 

end