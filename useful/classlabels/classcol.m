function [res] = classcol(labels)
%
% The function classcol generates a structure containing information on
% classes in matrix and vector formats. It also assigns colors to samples
% as a function of the class information. The colors are chosen among an
% internal predefined list of most distinguishable colors.
%
% Input argument :
% ================
% labels : vector containing samples' class information with any of the
% following formats : numeric, cell, table, string or character
%
% Output argument :
% =================
% res : structure containing unique labels, matrix/vector of samples'
% group information and matrix of colors for each sample
%
% Usage :
% =======
% [res] = classcol(labels);
%
% Related function :
% ==================
% distcol.m (generates distinguishable colors)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
%
% =========================================================================

% Converts cell, character and table formats to string 
if ischar(labels)==1 || iscell(labels)
    labels = string(labels);
    labels = strip(labels);
elseif istable(labels)==1
    labels = string(table2cell(labels));
end

% Identifies unique labels in the input
if isa(labels,'double')
    uni = unique(labels);
else
    uni = unique(labels,'stable');
end


% Sorts unique labels in alpha-numerical order
if isstring(labels) == 1 && sum(isnan(str2double(uni))) == 0
    uni = string(sort(str2double(uni)));
end

% Allocates space for classes and color matrices
classmat = zeros(length(labels),length(uni));
colmat = zeros(length(labels),3);

% Chooses distinguishable colors according to the number of unique labels
unicol = distcol(length(uni));

% If statement determining whether labels are numeric or string
% Fills classes and color matrices
if isnumeric(labels) == 1
    for i = 1 : length(uni)
        x = find(labels == uni(i));
        classmat(x,i) = 1;
        colmat(x,:) = repmat(unicol(i,:),[length(x) 1]);
    end
elseif isstring(labels) == 1
    for i = 1 : length(uni)
        x = find(strcmp(labels,uni(i)));
        classmat(x,i) = 1;
        colmat(x,:) = repmat(unicol(i,:),[length(x) 1]);
    end
end

% Converts groups matrix to vector format
classvec = classmat * (1:size(classmat,2))';

% Fills output structure
res.labels = uni;
res.labelscol = unicol;
res.classmat = classmat;
res.classvec = classvec;
res.colmat = colmat;

end
