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
% lme = fitlmematrix(x, Y, Z, G,...
%     'FixedEffectPredictors',...
%         {'ShapeParameter','ROI','Intercept'},...
%     'RandomEffectPredictors',...
%         {{'Intercept'}},...
%     'RandomEffectGroups',...
%         {'Subject'});


% CodeFolder = 'D:\github\AV-Attention-7T_code';
% FigureFolder = fullfile(CodeFolder, 'Figures');
% load(fullfile(FigureFolder, 'LMM_results.mat'), 'models');

load('LMM_results.mat', 'models');



%% diplay data and model
model = models(1);

figure(1)
colormap gray

subplot(1,4,1)
imagesc(model.Y)
set(gca, 'ytick', 1:numel(model.Y), 'yticklabel', model.Y_legend)
title('data')

subplot(1,4,2)
imagesc(model.X)
set(gca, 'xtick', 1:3, 'xticklabel', {'cst/lin'; 'ROI'; 'intercept'})
title('fixed effect design matrix')

subplot(1,4,3)
imagesc(model.Z)
set(gca, 'xtick', 1, 'xticklabel', {'intercept'})
title('random effect design matrix')

subplot(1,4,4)
imagesc(model.G)
set(gca, 'xtick', 1, 'xticklabel', {'subjects'})
title('random effect group')



%% print out results
model_of_interest = [1 4 3 2 5 8 9 10];

for i_model = model_of_interest %1:numel(models)
    
    model = models(i_model);
    
    % print out results
    fprintf('%s %i - %s', '%', i_model, model.name)
    fprintf(' - %s', model.ROIs)
    disp(model.lme)
    
    fprintf('\n\n\n')
    
end