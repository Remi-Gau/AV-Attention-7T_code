function display_lmm_results()% small script to print out the results of the LMM

clear
close all
clc

DO = 2; % 1 BOLD ; 2 MVPA

CodeFolder = 'D:\github\AV-Attention-7T_code';
FigureFolder = fullfile(CodeFolder, 'Figures');

% in case it is not running on remi's computer
if ~exist(FigureFolder, 'dir')
    FigureFolder = pwd;
end

NbSubj = 11;


% Models BOLD
% 1 - [A-fix]_{Att_A, Att_V} - A1 - PT
% 4 - [V-fix]_{Att_A, Att_V} - V1 - V2-3
% 3 - [V-fix]_{Att_A, Att_V} - A1 - PT
% 2 - [A-fix]_{Att_A, Att_V} - V1 - V2-3
% 5 - [AV - A]_{Att_A, Att_V} - A1 - PT
% 8 - [AV - V]_{Att_A, Att_V} - V1 - V2-3
% 9 - [Att_V - Att_A]_{A, V, AV} - A1 - PT
% 10 - [Att_V - Att_A]_{A, V, AV} - V1 - V2-3

% Models MVPA
% 1 - [AV VS A]_{att A, att V} - A1 - PT
% 4 - [AV VS V]_{att A, att V} - V1 - V2-3
% 5 - [Att_A VS Att_V]_{A, V, AV} - A1 - PT
% 6 - [Att_A VS Att_V]_{A, V, AV} - V1 - V2-3


% call :
% lme = fitlmematrix(x, Y, Z, G, 'FitMethod', 'REML',...
%     'FixedEffectPredictors',...
%     {'ROI1_cst', 'ROI1_lin', 'ROI2_cst', 'ROI2_lin'},...
%     'RandomEffectPredictors',...
%     {{'Intercept'}},...
%     'RandomEffectGroups',...
%     {'Subject'});

switch DO
    case 1
        load(fullfile(FigureFolder, 'LMM_BOLD_results.mat'), 'models');
        model_of_interest = [1 4 3 2 5 8 9 10];
    case 2
        load(fullfile(FigureFolder, 'LMM_MVPA_results.mat'), 'models');
        model_of_interest = [1 4:6];
end



%% diplay data and model
model = models(1);

figure(1)
colormap gray

% Data
subplot(1,4,1)
imagesc(model.Y)
set(gca, 'ytick', 1:numel(model.Y), 'yticklabel', model.Y_legend)
title('data')

% fixed effect design matrix
subplot(1,4,2)
imagesc(model.X)
set(gca, 'xtick', 1:4, 'xticklabel', {'ROI1 cst', 'ROI1 lin', 'ROI2 cst', 'ROI2 lin'})
title('fixed effect design matrix')

% random effect
% design matrix
subplot(1,4,3)
imagesc(model.Z)
set(gca, 'xtick', 1, 'xticklabel', {'intercept'})
title('random effect design matrix')

% group factor
subplot(1,4,4)
imagesc(model.G)
set(gca, 'xtick', 1, 'xticklabel', {'subjects'})
title('random effect group')



%% print out results
pattern = '%s\t F(%i,%i)= %f\t p = %f\n';
ToPermute = list_permutation();

for i_model = model_of_interest %1:numel(models)
    
    model = models(i_model);
    
    %print out results
    fprintf('%s %i - %s - %s \n', '%', i_model, model.name, model.ROIs)
    %     disp(model.lme)
    
    % display reults perm test and t-test for each s parameter for each ROI
    %     for i = 1:4
    %         compare_results(i, model, ToPermute);
    %     end
    %     fprintf('\n')
    
    % effect of either linear or constant in either ROIs
    c = [...
        1 0 0 0 ;...
        0 1 0 0 ;...
        0 0 1 0 ;...
        0 0 0 1 ];
    message = 'effect of either linear or constant in either ROI';
    PVAL = test_and_print(model, c, pattern, message);
    
    %  run LMM on just the CST or LIN from both ROIs if signiicant
    if PVAL<.05
        
        X = [...
            [ones(NbSubj,1) ; zeros(NbSubj,1)], ... % reg 1 : ROI1 cst
            [zeros(NbSubj,1) ; ones(NbSubj,1)]]; % reg 1 : ROI1 lin
        Z = repmat(ones(NbSubj,1), 2, 1);
        G = repmat([1:NbSubj]', 2, 1);
        
        % loop over cst and lin
        for i_s_param = 1:2
            
            switch i_s_param
                case 1
                    name_param = {'ROI1_cst', 'ROI2_cst'};
                    reg_of_interest = [1 3];
                case 2
                    name_param = {'ROI1_lin', 'ROI2_lin'};
                    reg_of_interest = [2 4];
                    
            end
            
            Y = model.Y(any(model.X(:,reg_of_interest), 2));
            
            submodel.lme = fitlmematrix(X, Y, Z, G, 'FitMethod', 'REML',...
                'FixedEffectPredictors',...
                name_param,...
                'RandomEffectPredictors',...
                {{'Intercept'}},...
                'RandomEffectGroups',...
                {'Subject'});
            
            c = [...
                1 0;...
                0 1];
            message = ['effect of ' name_param{1} ' and ' name_param{2} ' averaged across ROIs'];
            PVAL = test_and_print(submodel, c, pattern, message);
        end
        
    end
    
    fprintf('\n')
    
end

end


function [PVAL,F,DF1,DF2] = test_and_print(model, c, pattern, message)
[PVAL,F,DF1,DF2] = coefTest(model.lme, c);
fprintf(pattern, ...
    message, ...
    DF1, DF2, ...
    F, PVAL);
end

function ToPermute = list_permutation()
% create permutations for exact sign permutation test
for iSubj=1:11
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:}); clear sets
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
end

function Perms = create_null_distribution(ToPermute, Betas)
for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:)';
    Perms(iPerm,:) = mean(Betas.*repmat(tmp2,1,size(Betas,2))); %#ok<*AGROW>
end
end

function P = perm_test(betas, side, perms)
if strcmp(side,'left')
    P = sum(perms<mean(betas)) / numel(perms);
    
elseif strcmp(side,'right')
    P = sum(perms>mean(betas))/  numel(perms);
    
elseif strcmp(side,'both')
    P = sum(abs(perms)>abs(mean(betas))) / numel(perms) ;
end
end

function compare_results(i, model, ToPermute)

NbSubj = 11;
ROI_nb = [1 1 2 2];
side_idx = [1 2 1 2];
s_param = {'Cst', 'Lin', 'Cst', 'Lin'};

betas = model.Y(logical(model.X(:,i))); % get betas
side = model.test_side{side_idx(i)}; % get side for the test
perms = create_null_distribution(ToPermute, betas);

p_perm = perm_test(betas, side, perms);
[~, p_ttest] = ttest(betas, 0, 'tail', side);

% display the results of perm and t-test
fprintf('ROI %i - %s - p(perm) = %f - p(ttest) = %f\n', ...
    ROI_nb(i), s_param{i}, ...
    p_perm, p_ttest);

% OVERKILL: use LMM to do a t-test to make sure we get the same
% thing
lme = fitlmematrix(ones(NbSubj,1), betas, ones(NbSubj,1), [1:NbSubj]', 'FitMethod', 'REML',...
    'FixedEffectPredictors',...
    {'s_param'},...
    'RandomEffectPredictors',...
    {{'Intercept'}},...
    'RandomEffectGroups',...
    {'Subject'});
disp(lme.Coefficients)

end