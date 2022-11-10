function [mbsadecomp_res] = mbsadecomp(X, Y, adecomp_res)
%
% The function sadecomp stands for Synthetic ANOVA decomposition. It
% applies ANalysis Of VAriance decomposition on X according to a previous
% ANOVA decomposition as a reference.
%
% This function was first designed to be used when a rebalancing procedure
% is used before performing an ANOVA decomposition. It produces the mean
% effect matrices according to the ANOVA model along with the ANOVA
% matrices.
%
% Reference :
% ===========
% de Figueiredo, M., Giannoukos, S., Rudaz, S., Zenobi, R., & Boccard, J.
% (s. d.). Efficiently handling high-dimensional data from multifactorial
% designs with unequal group sizes using Rebalanced ASCA (RASCA).
% Journal of Chemometrics, n/a(n/a), e3401. https://doi.org/10.1002/cem.3401
%
% Input arguments :
% =================
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample
%     identified by integers (n x q), the synthetic design matrix
%
% adecomp_res : the ANOVA "model" built using the adecomp.m function
%
% Output arguments :
% ==================
% adecomp_res : structure containing results of the ANOVA decomposition
%   adecomp_res.Y: numerical array of the experimental design;
%   adecomp_res.Y: numerical array of the experimental data;
%   adecomp_res.Xm: the grand mean matrix;
%   adecomp_res.Xf: cell array of the pure effect matrices;
%   adecomp_res.Xe: numerical array of the pure error (residuals);
%   adecomp_res.ssq: sum of squares of the effects;
%   adecomp_res.ssqvarexp: sum of squares explained variation;
%   adecomp_res.Xfaug: cell array of the augmented effect matrices;
%   adecomp_res.effects: cell array of the effects (factor indices);
%   adecomp_res.Options: the options used for calculation;
%
% If the ANOVA decomposition uses the GLM methodology, also contains:
%   adecomp_res.Ef: the residuals without considering the effect f in the model;
%   adecomp_res.cod: the coding of the experimental design Y;
%   adecomp_res.B: the GLM parameters;
%
% Usage :
% =======
% [sadecomp_res] = sadecomp(X, Y, adecomp_res);
%
% Related functions :
% ===================
% adecomp.m (performs ANOVA decomposition)
% radecomp.m (performs rebalanced ANOVA decomposition)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
%
% =========================================================================

% Checks if columns of Y contain only consecutive positive integers
for i = 1 : size(Y,2)
    Y(:,i) = yformatconv(Y(:,i),'intvec');
end

ntab = length(X);
[n,~] = size(X{1}); % size of X
[~,q] = size(Y); % size of Y

% Loops over each table to perform a synthetic ANOVA decomposition
for t = 1 : ntab

    % Number of effects
    neff = length(adecomp_res(t).effects);

    % Allocates space for outputs
    Xf = cell(1,length(adecomp_res(t).Xf));
    Ef = cell(1,length(adecomp_res(t).Xf));
    ssq = zeros(1,neff);
    ssqvarexp = zeros(1,neff);

    % Retrieves options used for adecomp
    Options = adecomp_res(t).Options;

    switch Options.decomp

        case 'classical'

            % Retrieves the mean effect vector of each level, or combination of
            % levels, in each effect.
            Ysm = cell(1,neff);
            for i = 1 : neff

                Ys = adecomp_res(t).Y(:,adecomp_res(t).effects{i});

                Ysm{1,i}{1} = unique(Ys,'rows');

                for j = 1 : size(Ysm{1,i}{1},1)

                    [~,idx] = ismember(Ys,Ysm{1,i}{1}(j,:), 'rows');
                    Ysm{1,i}{2}(j,:) = mean(adecomp_res(t).Xf{i}(idx == 1, :),1);

                end

            end

            % Loops over the number of effects + the grand mean
            for i = 1 : neff + 1 % +1 for the grand mean

                if i == 1 % grand mean

                    Xm = repmat(adecomp_res(t).Xm(1,:),[size(Y,1),1]);

                else % effects

                    Xf{i-1} = zeros(size(X));
                    for k = 1 : size(Ysm{i-1}{1},1)
                        idx = ismember(Y(:, unique(adecomp_res(t).effects{i-1})), Ysm{i-1}{1}(k,:),'rows');
                        Xf{i-1}(idx,:) = repmat(Ysm{i-1}{2}(k,:),[sum(idx==1),1]);
                    end

                end

            end

            % Calculates the full residuals
            Xe = X - Xm - sum(cat(3,Xf{:}),3);

            % Calculates the sum of squares and explained variation
            for i = 1 : neff + 1

                if i < neff + 1 % effects

                    % Sum-of-squares
                    ssq(1,i) = sum(Xf{1,i}(:).^2);

                    % Explained variation
                    ssqvarexp(1,i) = ssq(1,i) / ( sum(X(:).^2) - sum(Xm(:).^2) ) * 100;

                else % residuals

                    % Sum-of-squares
                    ssq(1,i) = sum(Xe(:).^2);

                    % Explained variation
                    ssqvarexp(1,i) = ssq(1,i) / ( sum(X(:).^2) - sum(Xm(:).^2) ) * 100;

                end

            end

        case {'glm','rebalanced'}

            % Checks which coding is to be used to define the design matrices
            if strcmp(Options.coding,'sumcod') == 1

                % Sum coding of the input design matrix according to
                % Thiel et al. (2017)
                [cod, codempty] = sumcoding(Y, Options);

            elseif strcmp(Options.coding,'wecod') == 1

                % Weighted-effect coding of the input design matrix according
                % to Nieuwenhuis et al. (2017)
                [cod, codempty] = wecoding(Y, Options);

            end

            % Retrieves the B model parameters from adecomp_res
            B = adecomp_res(t).B;
            Bmat = cat(2,B{:});

            % Creates empty B useful for calculations
            bempty = cell(1,length(B));
            for i = 1 : length(B)
                bempty{i} = zeros(size(adecomp_res(t).X,2),size(cod{i},2));
            end

            % Calculates full model Residuals
            Xe = X{t} - (cat(2,cod{:}) * Bmat');

            % Calculation of the main factor effects :
            for i = 1 : length(cod) % all main effects and interactions

                % Calculates the effect matrices
                Xftmp = codempty; Xftmp{i} = cod{i};
                Bf = bempty; Bf{i} = B{i};

                if i == 1 % grand mean matrix

                    Xm = cat(2,Xftmp{:}) * cat(2,Bf{:})';

                else % effects

                    % Stores pure effect matrix
                    Xf{1,i-1} = cat(2,Xftmp{:}) * cat(2,Bf{:})';

                    % Calculates the residuals of the model with effect i
                    Xftmp = cod; Xftmp{i} = zeros(size(cod{i}));
                    Bf = B; Bf{i} = zeros(size(B{i}));
                    Ef{1,i-1} = X{t} - (cat(2,Xftmp{:}) * cat(2,Bf{:})');

                end

            end

            % Sum of squares and explained variance
            neff = length(cod)-1;
            den = sum(cat(3,Ef{1:neff}).^2,'all') - (neff-1) * sum(Xe(:).^2);
            for i = 1 : neff + 1

                if i < neff + 1 % effects

                    % Sum of squares
                    ssq(1,i) = sum(Ef{i}(:).^2) - sum(Xe(:).^2);

                    % Explained variance (row 1 does not sum to 100%, row 2 sums to 100%)
                    ssqvarexp(1,i) = ssq(1,i) / ( sum((X{t}-Xm).^2,'all') ) * 100;
                    ssqvarexp(2,i) = ssq(1,i) / ( den ) * 100;

                else % residuals

                    % Sum of squares
                    ssq(1,i) = sum(Xe(:).^2);

                    % Explained variance (row 1 does not sum to 100%, row 2 sums to 100%)
                    ssqvarexp(1,i) = sum(Xe(:).^2) / ( sum((X{t}-Xm).^2,'all') ) * 100;
                    ssqvarexp(2,i) = sum(Xe(:).^2) / ( den ) * 100;

                end

            end

    end

    % Calculates the augmented effect matrices Xfaug
    Xfaug = cell(1,length(Xf));
    for i = 1 : length(Xf)
        Xfaug{i} = Xf{i} + Xe;
    end

    % Prepares output structure
    sadecomp_res.Y = Y; % experimental design
    sadecomp_res.X = X{t}; % experimental data
    sadecomp_res.Xm = Xm; % grand mean
    sadecomp_res.Xf = Xf; % pure effects
    sadecomp_res.Xe = Xe; % residuals
    sadecomp_res.ssq = ssq; % sum of squares
    sadecomp_res.ssqvarexp = ssqvarexp; % explained variation
    sadecomp_res.Xfaug = Xfaug; % augmented effects

    if sum(strcmp(Options.decomp, {'glm','rebalanced'})) == 1
        sadecomp_res.Ef = Ef; % residuals without considering each effect f
        sadecomp_res.cod = cod; % glm coding of Y
        sadecomp_res.B = B; % glm parameters
    end

    sadecomp_res.effects = adecomp_res(t).effects; % effects
    sadecomp_res.Options = Options; % options used

    % Stores the results of the synthetic ANOVA decomposition
    mbsadecomp_res(t) = sadecomp_res;

end

end
