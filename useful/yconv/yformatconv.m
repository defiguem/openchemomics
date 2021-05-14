function [Yout] = yformatconv(Yin, ytarget)
% 
% The function yformatconv converts the input class assignment variable to 
% the format specified in the string ytarget. Supported formats are
% described below in input arguments.
% 
% Note that if input Yin is a vector, it should be in a column but more
% generally, samples are expected to be in the rows. Hence, the number of
% columns may vary if the format is binary matrix or not. Plus, trying to
% convert an input Yin with more than two classes to a two-class format
% like a binary vector or a vector of of ones and minus ones will throw an
% error.
% 
% Input arguments :
% =================
% Yin : class assignment variable with samples in rows, which can have
% one of the formats supported and defined by the input argument ytarget
% 
% ytarget : string defining of ouput format for class assignment variable
%   'binarymat' : binary matrix of 0s and 1s with samples in rows and classes in columns
%   'binaryvec' : binary vector of 0s and 1s with samples in rows and classes in one colum
%   'intvec' : column vector of consecutive integers (0s not allowed)
%   'pmones' : column vector of plus ones and minus ones
% 
% Output arguments :
% ==================
% Yout : output class assignment variable with conversion defined by ytarget
% 
% USAGE :
% =======
% ytarget = 'binarymat'; % or 'binaryvec' or 'intvec' or 'pmones'
% [Yout] = yformatconv(Yin,ytarget)
% 
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
% 
% Modifications:
% ==============
% 
% =========================================================================

%% Checks for NaNs

ydim = size(Yin); % dimensions of input Yin
nans = sum(isnan(Yin),2); % finds NaNs
Yin(nans > 0,:) = []; % removes these rows from Yin

%% Format checker for input Y

% Checks if the desired format for input Y conversion was specifided
if nargin < 2
    error('At least the targeted Y format must be specified (see yformatconv.m for supported formats)');
end

% If Y is a vector, checks if it is in column format
if isvector(Yin) == 1 && iscolumn(Yin) == 0
    Yin = Yin';
end

% Defines unique values in Y
uniquey = unique(Yin); 

% Checks if input Y is binary
isbin = isequal(uniquey,[0;1]);

% Checks if non-binary formats have zeros, which is prohibited 
if isbin == 0 && isvector(Yin) == 1 && any(Yin == 0)
    error('Class assignments in Y for non-binary input format cannot be zero');
end

if isbin == 1
    
    if isvector(Yin) == 1
        yin = 'binaryvec';
    else
        yin = 'binarymat';
    end
    
elseif isbin == 0 && isvector(Yin) == 1 && isequal(uniquey,[-1;1]) == 0
    
    yin = 'intvec';
    
elseif isbin == 0 && isvector(Yin) == 1 && isequal(uniquey,[-1;1]) == 1
    
    yin = 'pmones';
    
else
    
    Yout = Yin;
    
    % Replacing NaNs is important for consistency
    if sum(nans) ~= 0 % if there were any NaNs in the input, replaces them
        mask = NaN(ydim(1),size(Yout,2)); % create NaN mask based on inital positions
        mask(nans == 0,:) = Yout;
        Yout = mask; % replaces the NaNs
    end
    
    return;
    
end

% If the input format is the same as the output one, return input
if strcmp(yin,ytarget) == 1
    
    Yout = Yin;
    
    % Replacing NaNs is important for consistency
    if sum(nans) ~= 0 % if there were any NaNs in the input, replaces them
        mask = NaN(ydim(1),size(Yout,2)); % create NaN mask based on inital positions
        mask(nans == 0,:) = Yout;
        Yout = mask; % replaces the NaNs
    end
    
    return;
end

%% Format conversion according to format of input Y

switch ytarget
    
    case 'binarymat'
        
        if strcmp(yin,'binaryvec')
            
            Yout = [Yin, (Yin < 1)];
            
        elseif strcmp(yin,'intvec')
            
            Yout=zeros(size(Yin,1), max(Yin)); 
            Yout(sub2ind(size(Yout),1:numel(Yin),Yin'))=1;
            
        elseif strcmp(yin,'pmones')
            
            Yout = [(Yin == 1), (Yin == -1)];
            
        end
        
    case 'binaryvec'
        
        if strcmp(yin,'binarymat')
            if size(Yin,2) > 2
                error('Conversion to a binary vector cannot be performed if there are more than two-classes in input Y');
            else
                Yout = Yin(:,1);
            end
        elseif strcmp(yin,'intvec')
            if isequal(uniquey,[1;2]) == 0
                error('Conversion to a binary vector cannot be performed if there are more than two-classes in input Y');
            else
                Yout = Yin;
                Yout(Yout == 2) = 0;
            end
        elseif strcmp(yin,'pmones')
            Yout = double((Yin == 1));
        end
        
    case 'intvec'
        
        if strcmp(yin,'binarymat')
            Yout = Yin * (1:size(Yin,2))';
        elseif strcmp(yin,'binaryvec')
            Yout = Yin;
            Yout(Yout == 0) = 2;
        elseif strcmp(yin,'pmones')
            Yout = Yin;
            Yout(Yout == -1) = 2;
        end
        
    case 'pmones'
        
        if strcmp(yin,'binarymat')
            if size(Yin,2) > 2
                error('Conversion to a plus/minus one vector cannot be performed if there are more than two-classes in input Y');
            else
                Yout = Yin(:,1);
                Yout(Yout == 0) = -1;
            end
        elseif strcmp(yin,'binaryvec')
            Yout = Yin;
            Yout(Yout == 0) = -1;
        elseif strcmp(yin,'intvec')
            if isequal(uniquey,[1;2]) == 0
                error('Conversion to a binary vector cannot be performed if there are more than two-classes in input Y');
            else
                Yout = Yin;
                Yout(Yout == 2) = -1;
            end
        end
        
end

% Replacing NaNs is important for consistency
if sum(nans) ~= 0 % if there were any NaNs in the input, replaces them
    mask = NaN(ydim(1),size(Yout,2)); % create NaN mask based on inital positions
    mask(nans == 0,:) = Yout;
    Yout = mask; % replaces the NaNs
end

end