function [preproc, Options] = xpreproc(data, Options)
%
% The function xpreproc preprocesses data matrices with samples in the
% rows and variables in the columns. The list of the given preprocessing
% techniques available is given below and it is always possible to add new
% ones at will. For each preprocessing technique, the parameters used for
% the calculation are also stored in the output structure. This is
% important, especially when data was divided into calibration and
% validation sets. The parameters calculated with xpreproc for the
% calibration set are then used as inputs in the xpreprop function to
% preprocess the validation set. In that sense, each modification brought
% to the present function must be accompanied with the proper modification
% of the xpreprop function.
%
% Care must be taken in keeping the independence between the calibration
% and validation sets. In fact, some techniques such as the column mean
% centering use the parameters calculated with the calibration set to
% preprocess the validation set. Others, such as SNV work on the individual
% samples and each sample is preprocessed according to its own parameters.
% In the latter case, there is no difference between the calculations made
% with xpreproc or xpreprop.
%
% This function handles preprocessing of one or multiple blocks. Usage must
% be adapted accordingly. For each table, one single or several
% pretreatments can be performed in a row. See usage below for more
% information.
%
% References :
% ============
% Roger, J.-M., Boulet, J.-C., Zeaiter, M., & Rutledge, D. N. (2020).
% 3.01-Pre-processing Methods. In S. Brown, R. Tauler, & B. Walczak (Éds.),
% Comprehensive Chemometrics (Second Edition) (p. 1-75). Elsevier.
% https://doi.org/10.1016/B978-0-12-409547-2.14878-4
%
% Input arguments :
% =================
% data : matrix of samples in the rows and variables in the columns;
%
% Options.prepro.X.type : array of cells containing in each cell the chain
% of pretreatments to apply to each data block;
%
% Options.prepro.X.para : array of cells containing in each cell the
% parameters associated to each pretreatment in Options.prepro.X.type;
%
% Output arguments :
% ==================
% preproc : structure containing the preprocessed data in preproc.data
% along with sub-structures for each of the preprocessing techniques used
% containing the parameters used for calculation.
%
% Usage :
% =======
% Options.prepro.X.type = {'cmeancenter'};
% Options.prepro.X.para = repelem({''},length(Options.prepro.X.type));
% [preprocx] = xpreproc(data, Options);
%
% or
%
% Options.prepro.X.type = {'snv';'cmeancenter'};
% Options.prepro.X.para = repelem({''},length(Options.prepro.X.type));
% [preprocx] = xpreproc(data, Options);
%
% or
%
% Options.prepro.X.type = {{'autoscale'};{'cmeancenter'}};
% Options.prepro.X.para = repelem({''},length(Options.prepro.X.type));
% [preprocx] = xpreproc(data, Options); % data is a cell array of 2 tables
%
% Related functions :
% ===================
% xpreprop.m (does preprocessing on X testing data)
% ypreproc.m (does preprocessing on Y training data)
% yunprepro.m (does unpreprocessing on Y preprocessed data)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications :
% ===============
%
% =========================================================================

%% Fail-safe section

% Checks if an array of cells with several tables was provided
if iscell(data) == 1
    ntab = length(data); % number of tables in data input
else
    data = {data}; % place array in cell
    ntab = 1;
end

% Checks if preprocessing methods were chosen for X, otherwise does nothing
if isfield(Options,'prepro') == 0 || isfield(Options.prepro,'X') == 0 || isfield(Options.prepro.X,'type') == 0
    Options.prepro.X.type = repelem({'nothing'},ntab)';
end

% Checks if preprocessing parameters were provided when needed
if isfield(Options.prepro.X,'para') == 0 
    for i = 1 : ntab
        Options.prepro.X.para{1,i} = repelem({''},length(Options.prepro.X.type))';
    end
end

% Checks if preprocessing methods were chosen for X, otherwise does nothing
if isfield(Options.prepro.X,'type') == 1 && ntab ~= length(Options.prepro.X.type)
    error('Preprocessing methods to be applied must be defined for each of the input tables');
end

if ischar(Options.prepro.X.type) == 1 % Cell format needed in the for loop
    Options.prepro.X.type = cellstr(Options.prepro.X.type);
end

%% Main section

% Loops over the block in the input argument data
for t = 1 : ntab
    
    if ischar(Options.prepro.X.type{t}) == 1 % Cell format needed in the for loop
        Options.prepro.X.type{t} = cellstr(Options.prepro.X.type{t});
    end
    
    % Stores the preprocessing steps to be used in the order of application
    preproc(t).preprosteps = Options.prepro.X.type{t};
    
    % Allocates space for storage of preprocessing methods parameters
    preproc(t).prepropara = cell( length(preproc(t).preprosteps), 1 );
    
    % Starts the preprocessing loop with as many iterations as preprocessing
    % techniques specified in Options.prepro.X.type{t}
    for i = 1 : length(Options.prepro.X.type{t})
        
        % Determines the dimensions of the data matrix
        % Needed inside the loop in case some preprocessing techniques modifiy the data size
        [xrow,xcol] = size(data{t});
        
        switch Options.prepro.X.type{t}{i} % Selects one preprocessing technique at a time
            
            case 'nothing'
                % No preprocessing performed
                
            case 'cmeancenter' % Column mean centering (no parameters to optimize)
                cmdata = mean(data{t},1,'omitnan'); % Column mean calculation
                data{t} = data{t} - ( ones(xrow,1) * cmdata ); % Substracts from each element in a given column the associated column mean
                
                % Parameter storage
                preproc(t).prepropara{i}.cmeancenter.cmdata = cmdata;
                
            case 'rmeancenter' % Row mean centering (no parameters to optimize)
                rmdata = mean(data{t},2,'omitnan'); % Row mean calculation
                data{t} = data{t} - ( rmdata * ones(1,xcol) ); % Substracts from each element in a given row the associated row mean
                
            case 'autoscale' % Column autoscaling / variance scaling / column standardization / z-transformation (no parameters to optimize)
                cmdata = mean(data{t},1,'omitnan'); % Column mean calculation
                cstddata = std(data{t},0,1,'omitnan'); % Column standard deviation calculation
                data{t} = (data{t} - ( ones(xrow,1) * cmdata )) ./ ( ones(xrow,1) * cstddata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation
                
                % Parameters storage
                preproc(t).prepropara{i}.autoscale.cmdata = cmdata;
                preproc(t).prepropara{i}.autoscale.cstddata = cstddata;
                
                % Autoscaling in the row dimension is called SNV (Standard Normal Variates)
                
            case 'robautoscale' 
                cmdata = median(data{t},1,'omitnan'); % Column median calculation
                ciqrdata = prctile(data{t}, 75) - prctile(data{t}, 25) ; % interquartile range of each row
                
                data{t} = (data{t} - ( ones(xrow,1) * cmdata )) ./ ( ones(xrow,1) * ciqrdata );
                
                % Parameters storage
                preproc(t).prepropara{i}.robautoscale.cmdata = cmdata;
                preproc(t).prepropara{i}.robautoscale.ciqrdata = ciqrdata;
                
            case 'stdscale' % Column autoscaling / variance scaling / column standardization / z-transformation (no parameters to optimize)
                
                    cstddata = std(data{t},0,1,'omitnan'); % Column standard deviation calculation
                    cstddata(cstddata <= 1e8) = 1;
                    
                    data{t} = data{t} ./ ( ones(xrow,1) * cstddata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation
                    
                    % Parameters storage
                    preproc(t).prepropara{i}.stdscale.cstddata = cstddata;

            case 'varscale' % Column autoscaling / variance scaling / column standardization / z-transformation (no parameters to optimize)
                
                    cvardata = var(data{t},0,1,'omitnan'); % Column standard deviation calculation
                    cvardata(cvardata <= 1e8) = 1;
                    
                    data{t} = data{t} ./ ( ones(xrow,1) * cvardata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation
                    
                    % Parameters storage
                    preproc(t).prepropara{i}.varscale.cvardata = cvardata;
                    
            case 'robscale'
                
                ciqrdata = prctile(data{t}, 75) - prctile(data{t}, 25) ; % interquartile range of each row
                %ciqrdata = mad(data{t}) ;
                
                data{t} = data{t} ./ ( ones(xrow,1) * ciqrdata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation

                % Parameters storage
                preproc(t).prepropara{i}.robscale.ciqr = ciqrdata;
                
            case 'cpareto' % Column Pareto scaling (no parameters to optimize)
                cmdata = mean(data{t},1); % Column mean calculation
                csqstddata = sqrt(std(data{t},0,1)); % Column square root of standard deviation calculation
                data{t} = ( data{t} - ( ones(xrow,1) * cmdata ) ) ./ ( ones(xrow,1) * csqstddata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation square root
                
                % Parameters storage
                preproc(t).prepropara{i}.cpareto.cmdata = cmdata;
                preproc(t).prepropara{i}.cpareto.csqstddata = csqstddata;
                
            case 'paretoscale' % Column Pareto scaling (no parameters to optimize)
                csqstddata = sqrt(std(data{t},0,1)); % Column square root of standard deviation calculation
                data{t} = data{t} ./ ( ones(xrow,1) * csqstddata ); % Substracts from each element in a given column the associated column mean and divides it by the associated column standard deviation square root
                
                % Parameters storage
                preproc(t).prepropara{i}.cpareto.csqstddata = csqstddata;
                
            case 'rpareto' % Row Pareto scaling (no parameters to optimize)
                rmdata = mean(data{t},2); % Row mean calculation
                rsqstddata = sqrt( std(data{t},0,2) ); % Row standard deviation square root calculation
                data{t} = (data{t} - ( rmdata * ones(1,xcol) ) ) ./ ( rsqstddata * ones(1,xcol) ); % Substracts from each element in a given row the associated row mean and divides it by the associated row standard deviation square root
                
            case 'cmaxvalnorm' % Column Max Normalization (no parameters to optimize)
                cmaxval = max(data{t},[],1); % max value for each column
                data{t} = data{t} ./ (ones(xrow,1) * cmaxval); % max normalization
                
                % Parameters storage
                preproc(t).prepropara{i}.cmaxvalnorm.cmaxval = cmaxval;
                
            case 'rmaxvalnorm' % Row Max Normalization (no parameters to optimize)
                rmaxval = max(data{t},[],2); % max value for each row
                data{t} = data{t} ./ (rmaxval * ones(1,xcol)); % max normalization
                
            case 'cnorm' % Column normalization (no parameters to optimize)
                cnormval = sqrt(sum(data{t}.^2,1,'omitnan')); % norm value for each column
                data{t} = data{t} ./ (ones(xrow,1) * cnormval); % normalization
                
                % Parameters storage
                preproc(t).prepropara{i}.cnorm.cnormval = cnormval;
                
            case 'rnorm' % Row normalization (no parameters to optimize)
                normval = sqrt(sum(data{t}.^2,2,'omitnan')); % norm value for each row
                data{t} = data{t} ./ (normval * ones(1,xcol)); % normalization
                
            case 'cminmaxnorm' % Column Min-Max Range Normalization (no parameters to optimize)
                cminmax = max(data{t},[],1,'omitnan') - min(data{t},[],1,'omitnan'); % Min-Max range

                if any(cminmax == 0)
                    %idx = (cminmax == 0);
                    %cminmax(idx) = max(data{t}(:,idx),[],1,'omitnan');
                    cminmax = ones(1,length(cminmax));
                end

                data{t} = data{t} ./ (ones(xrow,1) * cminmax); % Min-Max range norm

                % Parameters storage
                preproc(t).prepropara{i}.cminmaxnorm.cminmax = cminmax;
                
            case 'rminmaxnorm' % Row Min-Max Range Normalization (no parameters to optimize)
                minmax = max(data{t},[],2) - min(data{t},[],2); % Min-Max range
                data{t} = data{t} ./ (minmax * ones(1,xcol)); % Min-Max range norm
                
            case 'cminmaxcorr' % Column Min-Max correction (no parameters to optimize)
                if sum(isfield(Options.prepro,{'desmin','desmax'})) == 2
                    desmin = Options.prepro.desmin; % desired min value
                    desmax = Options.prepro.desmax; % desired max value
                else
                    desmin = 0;
                    desmax = 1;
                    Options.prepro.desmin = desmin; % desired min value
                    Options.prepro.desmax = desmax; % desired max value
                end
                
                data{t} = (data{t} - ( ones(xrow,1) * min(data{t},[],1)) );
                data{t} = data{t} ./ ( ones(xrow,1) * max(data{t},[],1) ); % normalized first to range 0-1
                data{t} = (data{t} .* ( ones(xrow,1) * (desmax - desmin) ) ) + ( ones(xrow,1) * desmin ); % scaled to desired min/max range
                
                % Parameters storage
                preproc(t).prepropara{i}.cminmaxcorr.desmin = desmin;
                preproc(t).prepropara{i}.cminmaxcorr.desmax = desmax;
                
            case 'rminmaxcorr' % Row Min-Max correction (no parameters to optimize)
                if sum(isfield(Options.prepro,{'desmin','desmax'})) == 2
                    desmin = Options.prepro.desmin; % desired min value
                    desmax = Options.prepro.desmax; % desired max value
                else
                    desmin = 0;
                    desmax = 1;
                    Options.prepro.desmin = desmin; % desired min value
                    Options.prepro.desmax = desmax; % desired max value
                end
                
                data{t} = (data{t} - (min(data{t},[],2) * ones(1,xcol)) );
                data{t} = data{t} ./ ( max(data{t},[],2) * ones(1,xcol) ); % normalized first to range 0-1
                data{t} = (data{t} .* ( (desmax - desmin) * ones(1,xcol)) ) + ( desmin * ones(1,xcol) ); % scaled to desired min/max range
                
                % Parameters storage
                preproc(t).prepropara{i}.rminmaxcorr.desmin = desmin;
                preproc(t).prepropara{i}.rminmaxcorr.desmax = desmax;
                
            case 'log10' % logarithmic transformation (no parameters to optimize)
                data{t}(data{t} ~= 0) = real( log10(data{t}(data{t} ~= 0)) );

            case 'ln' % logarithmic transformation (no parameters to optimize)
                data{t}(data{t} ~= 0) = real( log(data{t}(data{t} ~= 0)) );

            case 'simpleratios'
                indices = nchoosek(1:size(data{t},2),2);
                data{t} = data{t}(:,indices(:,1)) ./ data{t}(:,indices(:,2));

                % Parameters storage
                preproc(t).prepropara{i}.simpleratios.indices = indices;

            case 'doubleratios'
                indices = nchoosek(1:size(data{t},2),2);
                data{t} = data{t}(:,indices(:,1)) ./ data{t}(:,indices(:,2));

                % Parameters storage
                preproc(t).prepropara{i}.doubleratios.indices = indices;
                
            case 'snv' % Standard Normal Variates (no parameters to choose)
                rmdata = mean(data{t},2); % Row mean calculation
                rstddata = std(data{t},0,2); % Row standard deviation calculation
                data{t} = ( data{t} - rmdata * ones(1,xcol) ) ./ ( rstddata * ones(1,xcol) ); % Substracts from each element in a given row the associated row mean and divides it by the associated row standard deviation
                
            case 'rnv' % Robust Normal Variates (Roger et al., 2020) (may be needed to optimize the percentile value to use)
                if isfield(Options.prepro,'rnvprct')==0
                    Options.prepro.rnvprct = 10;
                end
                
                prct = Options.prepro.rnvprct;
                
                xprct = prctile(data{t}, prct, 2); % percentile of each row
                data{t} = ( data{t} - xprct * ones(1,xcol) ); % row robust center
                
                xiqr = abs( prctile(data{t}, prct, 2) - prctile(data{t}, (100 - prct), 2) ); % interquartile range of each row
                data{t} = data{t} ./ (xiqr * ones(1,xcol)); % row robust standardization
                
                % Parameters storage
                preproc(t).prepropara{i}.rnv.prct = prct;
                
%             case 'vsn' % Variable Sorting for Normalization : Rabatel, G., Marini, F., Walczak, B., & Roger, J.-M. (2020). VSN?: Variable sorting for normalization. Journal of Chemometrics, 34(2), e3164. https://doi.org/10.1002/cem.3164
%                 if isfield(Options.prepro.X.para{t},'vsnnpar') == 0
%                     Options.prepro.X.para{t}.vsnnpar = 2;
%                 end
%                 
%                 if isfield(Options.prepro.X.para{t},'vsntol') == 0
%                     Options.prepro.X.para{t}.vsntol = 0.0067;
%                 end
%                 
%                 [data{t},res] = vsn(data{t},struct('nparameters',Options.prepro.X.para{t}.vsnnpar,'tolerance',Options.prepro.X.para{t}.vsntol));
%                 
%                 % Parameters storage
%                 preproc(t).prepropara{i}.vsn.res = res;
                
            
            case 'msc' % Multiplicative Scatter Correction (can give reference data{t} in Options.prepro.mscref, otherwise uses mean of matrix)
                if isfield(Options.prepro,'mscref') == 0
                    dataRef = mean(data{t},1);
                end
                
                for j = 1 : xrow
                    b = [ones(xcol,1),dataRef'] \ data{t}(j,:)';
                    data{t}(j,:) = (1/b(2)) .* ( data{t}(j,:) - b(1) * ones(1,xcol) );
                end
                
                % Parameters storage
                preproc(t).prepropara{i}.msc.ref = dataRef;
                
            case 'emsc' % Extended Multiplicative Scatter Correction (can give reference data{t} in Options.prepro.mscref, otherwise uses mean of matrix)
                if isfield(Options.prepro,'emscref') == 0
                    dataRef=mean(data{t},1);
                end
                
                lambda = (1:xcol)';
                lambda2 = lambda.^2;
                
                for j = 1 : xrow
                    b = [ones(xcol,1),dataRef',lambda,lambda2] \ data{t}(j,:)';
                    data{t}(j,:) = (1/b(2)) .* ( data{t}(j,:) - b(1) * ones(1,xcol) - b(3) * lambda' - b(4) * (lambda2') );
                end
                
                % Parameters storage
                preproc(t).prepropara{i}.emsc.ref = dataRef;
                
            case 'pqn' % Probabilistic Quotient Normalization
                % Additive effects, if present, must be corrected before applying PQN
                if isfield(Options.prepro,'pqnref') == 0
                    pqnref = median(data{t},1);
                end
                
                ratio = abs( data{t} ./ ( ones(xrow,1) * pqnref ));
                medratio = median(ratio,2);
                data{t} = data{t} ./ ( medratio * ones(1,xcol) );
                
                % Parameters storage
                preproc(t).prepropara{i}.pqn.ref = pqnref;
                
            case 'detrend' % No parameters to choose (can give reference data in Options.prepro.pqnref, otherwise uses mean of matrix), (New code Roger et al., 2020)
                %             % Old code
                %             lambda=(1:xcol)';
                %             lambda2=lambda.^2;
                %             for j=1:xrow
                %                 b=[ones(xcol,1),lambda,lambda2]\data(j,:)';
                %                 data(j,:)=(data(j,:)-b(1)*ones(1,xcol)-b(2)*lambda'-b(3)*(lambda2'));
                %             end
                if isfield(Options.prepro,'dpoly') == 0
                    Options.prepro.dpoly = 2;
                    dpoly = Options.prepro.dpoly;
                else
                    dpoly = Options.prepro.dpoly;
                end
                
                % Wavelengths vector
                wl = (1:xcol)';
                
                % Polynomial basis for othogonal projection
                lambda = ones(xcol,1);
                for j = 2 : dpoly
                    lambda = [lambda,wl.^j];
                end
                
                % Orthogonal projection of the data into the basis
                data{t} = data{t} - data{t} * lambda * pinv(lambda'*lambda) * lambda';
                
                % Parameters storage
                preproc(t).prepropara{i}.detrend.poly = dpoly;
                
            case 'spline' % Spline baseline correction (Roger et al., 2020)
                if isfield(Options.prepro,'splineap') == 0
                    error('At least anchor-points indices must be provided in Options.prepro.ap');
                else
                    ap = Options.prepro.splineap;
                end
                
                % Baselines computation for each sample
                splinetrend = zeros(xrow,xcol);
                for j = 1 : xrow
                    splinetrend(j,:) = spline(ap,data{t}(j,ap),[1:xcol]');
                end
                
                % Corrects all samples
                data{t} = data{t} - splinetrend;
                
                % Parameters storage
                preproc(t).prepropara{i}.spline.ap = ap;
                
            case 'als' % Asymmetric Least Sqaures (Roger et al., 2020)
                if isfield(Options.prepro,'alslambda') == 0
                    Options.prepro.alslambda = 1;
                    lambda = Options.prepro.alslambda;
                else
                    lambda = Options.prepro.alslambda;
                end
                
                if isfield(Options.prepro,'alsq') == 0
                    Options.prepro.alsq = 0.1;
                    q = Options.prepro.alsq;
                else
                    q = Options.prepro.alsq;
                end
                
                bl = zeros(xrow,xcol);
                for j = 1 : xrow
                    x = data{t}(j,:)';
                    delta = diff(speye(xcol),2);
                    w = ones(xcol,1);
                    w0 = 2;
                    for k = 1 : 1000 % arbitrary max number of iterations
                        % Solves minimization
                        W = spdiags(w,0,xcol,xcol);
                        C = chol(W + lambda * (delta' * delta));
                        z = C \ (C' \ (w .* x));
                        % Updates the weights
                        w = q * (x > z) + (1 - q) * (x < z);
                        % Convergence criterion
                        fit = sum( abs(w0 - w) );
                        if fit < 1e-10
                            break;
                        end
                        w0 = w;
                    end
                    bl(j,:) = z';
                end
                
                data{t} = data{t} - bl;
                
                % Parameters storage
                preproc(t).prepropara{i}.als.lambda = lambda;
                preproc(t).prepropara{i}.als.q = q;
                
            case 'movav' % Smoothing by moving average (Roger et al., 2020)
                if isfield(Options.prepro,'movavwin') == 0
                    Options.prepro.X.para{t}.movavwin = 11;
                    win = Options.prepro.X.para{t}.movavwin;
                else
                    win = Options.prepro.X.para{t}.movavwin;
                end
                
                w = floor(win/2);
                x = [data{t}(:,w:-1:1), data{t}, data{t}(:,end-w+1:end)];
                for j = w+1 : w+xcol
                    data{t}(:,j-w) = (sum(x(:,j-w:j+w)'))'./(2*w+1);
                end
                
                % Parameters storage
                preproc(t).prepropara{i}.movav.win = win;
                
            case 'savgol' % Need to choose window size plus degree for derivative and polynomial
                % Warning : a 2nd derivative needs at 2nd degree polynomial at least
                % and the polynomial has to be smaller than the odd window size
                % from at least two units
                
                
                poly = Options.prepro.X.para{t}.sgpoly;
                win = Options.prepro.X.para{t}.sgwin;
                deriv = Options.prepro.X.para{t}.sgderiv;
                
                [~,g] = sgolay(poly,win);
                dt = 1;
                
                tails = [data{t}(:,1 : win); data{t}(:,end-win+1 : end)];

%                % Iterative code over samples
%                for j = 1 : xrow
%                    data{t}(j,:) = conv(data{t}(j,:)',factorial(deriv)/(-dt)^deriv*g(:,deriv+1),'same')';
%                end
                
                data{t} = conv2(data{t}', factorial(deriv)/(-dt)^deriv * g(:,deriv+1), 'same')';
                
                % Handling preprocessing of endpoints/tails
                w = floor(win/2);
                x = repmat((-w:w)',1,1+poly).^repmat((0:poly), win,1);
                a = ((x'*x)\x')';
                a = tails * a;  
                for k = 1 : deriv
                    a = a(:,2 : poly+2-k) * diag(1 : poly+1-k); % or its d'th derivative
                end
                data{t}(:,1 : w+1) = a(1:size(data{t},1),:) * x(1:w+1,1:1+poly-deriv)'; 
                data{t}(:,end-w:end) = a(size(data{t},1)+(1:size(data{t},1)),:) * x(w+1:win,1:1+poly-deriv)';
                
                % Parameters storage
                preproc(t).prepropara{i}.savgol.poly = poly;
                preproc(t).prepropara{i}.savgol.win = win;
                preproc(t).prepropara{i}.savgol.deriv = deriv;
                
            case 'pct' % Principal Components Transform filtering (Roger et al., 2020)
                if isfield(Options.prepro,'pctpcs') == 0
                    error('Principal components to retain must be defined in a vector in Options.prepro.pctpcs');
                elseif Options.prepro.pctpcs >= min([xrow,xcol])
                    error('The number of principal components to filter is greater or equal to the total number of principal components available');
                end
                
                pcs = Options.prepro.pctpcs; % pcs to select
                
                [u,s,v] = svd(data{t}, 'econ');
                v = v(:,pcs); % retains components defined in pcs
                
                data{t} = data{t} * (v * v'); % projects data{t} into new basis
                
                % Parameters storage
                preproc(t).prepropara{i}.pct.loadings = v; % basis to project test samples
                
%             case 'ica' % Independent Components Analysis filtering (Roger et al., 2020)
%                 if isfield(Options.prepro,'icaics') == 0
%                     error('Independent components to retain must be defined in a vector in Options.prepro.icaics');
%                 elseif Options.prepro.icaics >= min([xrow,xcol])
%                     error('The number of independent components to filter is greater or equal to the total number of independent components available');
%                 end
%                 
%                 ics = Options.prepro.icaics; % ics to retain
%                 sg = jadeR(data{t}, max(ics)) * data{t}; % pure signals extracted
%                 sc = data{t} * sg' * pinv(sg * sg'); % calculates scores
%                 
%                 data{t} = sc * sg; % calculates filtered data{t}
%                 
%                 % Parameters storage
%                 preproc(t).prepropara{i}.ica.signals = sg; % stores pure signals
                
            case 'epo' % External Parameter Orthogonalization (Roger et al., 2020)
                if isfield(Options.prepro,'yref') == 0
                    error('Class association of samples must be defined in Options.prepro.yref as a consecutive integer vector');
                else
                    yref = Options.prepro.yref;
                end
                
                if isfield(Options.prepro,'epopcs') == 0
                    error('The number of principal components to retain must be defined when usin External Parameter Orthogonalization');
                end
                
                D = [];
                
                uni = unique(yref);
                for i = 1 : length(uni)
                    x = data{t}(yref == uno(i),:);
                    m = mean(x);
                    D = [D ; ( x - repmat(m,size(x,1),1) )];
                end
                
                [u,s,P] = svd(D, 'econ');
                P = P(:,Options.prepro.epopcs);
                
                data{t} = data{t} - (data{t} * P * P');
                
                % Parameters storage
                preproc(t).prepropara{i}.epo.loadings = P;
                
            case 'tabnorm' % normalization of table to unit variance
                tabnorm = sqrt( sum( sum( data{t}.^2 ) ) ); % frobenius norm
                data{t} = data{t} ./ tabnorm;
                
                % Parameters storage
                preproc(t).prepropara{i}.tabnorm = tabnorm;
                
            otherwise
                error('Preprocessing technique ''%s'' does not exist or is not supported', Options.prepro.X.type{t}{i});
                
        end
        
    end
    
    % The data is preprocessed accross the loop iterations and the final form
    % is stored in the output structure only in the end. However, the
    % parameters calculated throughout the loop iterations are directly stored.
    preproc(t).data = data{t};
    
end

end