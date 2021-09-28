function [T, P, ord] = pcsignord(T, P, Pref)
%
% Applying PCA to different subsets of the same dataset may cause 
% inversions in the order of the principal components extracted. Moreover, 
% the sign of the principal components being arbitrary in PCA, matching 
% components in two models may be of different signs, where one component 
% is the reflection of the other. Correcting for the order and the sign 
% ambiguity between two PCA model can be done through the calculation of a 
% correlation matrix between the PCs in the models.The re-ordering of the 
% PCs is based on maximum absolute correlation between the PCs 
% (independently of their sign), then the sign of theses correlations can 
% be used to determine if a reflection is present. If so, then the signs of 
% one of the PCs (scores and loadings) is multiplied by -1.
%
% References :
% ============
% Babamoradi, H., van den Berg, F., & Rinnan, Å. (2013). Bootstrap based 
% confidence limits in principal component analysis?A case study. 
% Chemometrics and Intelligent Laboratory Systems, 120, 97?105. 
% https://doi.org/10.1016/j.chemolab.2012.10.007
%
% Peres-Neto, P. R., Jackson, D. A., & Somers, K. M. (2003). Giving 
% Meaningful Interpretation to Ordination Axes?: Assessing Loading 
% Significance in Principal Component Analysis. Ecology, 84(9), 2347?2363. 
% https://doi.org/10.1890/00-0634
%
% Input arguments :
% =================
% T : scores of the PCA model to correct;
% P : loadings of the PCA model to correct;
% Pref : reference loadings towards which correction is made;
% 
% Output arguments :
% ==================
% T : corrected scores of the PCA model;
% P : corrected loadings of the PCA model;
% ord : order of the PCs to check  if reordering was made;
% 
% Usage :
% =======
% [T, P, ord] = pcsignord(T, P, Pref);
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
% 
% Modifications:
% ==============
%
% =========================================================================

% Corrects for ordering and sign ambiguities
R = corr(P,Pref); % calculates correlations between columns of P and Pref
Pcorr = P; Tcorr = T; % duplicates initial scores and loadings
ord = [];

for i = 1 : size(R,2) % loops over the columns
    
    [~,idx] = sort(abs(R(i,:)),'descend'); % find column in Pref, which has greatest correlation with of P(:,i)
    C = setdiff(idx,ord,'legacy');
    ord(i) = C(1);
    
    Pcorr(:,i) = P(:,ord(i)); % orders the loadings
    Tcorr(:,i) = T(:,ord(i)); % orders the scores
    if R(i,ord(i)) < 0 % corrects for sign ambiguity if correlation is negative
        Pcorr(:,i) = Pcorr(:,i) .* -1;
        Tcorr(:,i) = Tcorr(:,i) .* -1;
    end
end

T = Tcorr; 
P = Pcorr;
clear Tcorr Pcorr

end