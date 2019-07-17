% small script to print out the results of the LMM

clear
close all
clc

% Models
% 1 - [A-fix]_{Att_A, Att_V} - A1 - PT
% 4 - [V-fix]_{Att_A, Att_V} - V1 - V2-3
% 3 - [V-fix]_{Att_A, Att_V} - A1 - PT
% 2 - [A-fix]_{Att_A, Att_V} - V1 - V2-3
% 5 - [AV - A]_{Att_A, Att_V} - A1 - PT
% 8 - [AV - V]_{Att_A, Att_V} - V1 - V2-3
% 9 - [Att_V - Att_A]_{A, V, AV} - A1 - PT
% 10 - [Att_V - Att_A]_{A, V, AV} - V1 - V2-3

% call :
% lme = fitlmematrix(x, Y, Z, G, 'FitMethod', 'REML',...
%     'FixedEffectPredictors',...
%     {'ROI1_cst', 'ROI1_lin', 'ROI2_cst', 'ROI2_lin'},...
%     'RandomEffectPredictors',...
%     {{'Intercept'}},...
%     'RandomEffectGroups',...
%     {'Subject'});


% CodeFolder = 'D:\github\AV-Attention-7T_code';
% FigureFolder = fullfile(CodeFolder, 'Figures');
% load(fullfile(FigureFolder, 'LMM_results.mat'), 'models');

load('LMM_results.mat', 'models');



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
model_of_interest = [1 4 3 2 5 8 9 10];
pattern = '%s\t F(%i,%i)= %f\t p = %f\n';

for i_model = model_of_interest %1:numel(models)
    
    model = models(i_model);
    
    % print out results
    fprintf('%s %i - %s - %s', '%', i_model, model.name, model.ROIs)
    disp(model.lme)
    
    % results from some specific contrasts
    
    % effect of CST
    fprintf('%s %i - %s - %s\n', '%', i_model, model.name, model.ROIs)
    c = [1 0 1 0];
    [PVAL,F,DF1,DF2] = coefTest(model.lme, c);
    fprintf(pattern, ...
        'effect of CST', ...
        DF1, DF2, ...
        F, PVAL);
    
    % effect of LIN
    c = [0 1 0 1];
    [PVAL,F,DF1,DF2] = coefTest(model.lme, c);
    fprintf(pattern, ...
        'effect of LIN', ...
        DF1, DF2, ...
        F, PVAL);
    
    % effect of either
    c = [...
        1 0 1 0 ;...
        0 1 0 1 ];
    [PVAL,F,DF1,DF2] = coefTest(model.lme, c);
    fprintf(pattern, ...
        'effect of either CST/LIN', ...
        DF1, DF2, ...
        F, PVAL);
    
    
    fprintf('\n\n\n')
end