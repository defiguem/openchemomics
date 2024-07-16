function [mkr] = distmkr(n)
%
% Personal list of distinguishable markers sorted in preference order
%
% Input argument :
% ================
% n : number of distinguishable markers to extract
%
% Output argument :
% =================
% mkr : cell array with as markers as the input argument n 
%
% Usage :
% =======
% [mkr] = distmkr(n);
%
% Related function :
% ==================
% classcol.m 
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
%
% =========================================================================

mkrs = {'o'; ... % Circle
    'square'; ... % Square
    '^'; ... % Upward-pointing triangle
    'diamond'; ... % Diamond
    'pentagram'; ... % Five-pointed star (pentagram)
    'hexagram'; ... % Six-pointed star (hexagram)
    'v'; ... % Downward-pointing triangle
    '>'; ... % Right-pointing triangle
    '<'}; ... % Left-pointing triangle

% If n is greater than the type of markers available, marker assignment
% starts again from the beginning of the list
reps = ceil(n/length(mkrs));
idx = repmat(1:length(mkrs),[1,reps])';
mkr = mkrs(idx(1:n));

end