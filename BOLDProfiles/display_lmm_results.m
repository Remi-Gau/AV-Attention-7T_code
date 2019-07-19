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

for i_model = model_of_interest %1:numel(models)
    
    model = models(i_model);
    
    % print out results
    %     fprintf('%s %i - %s - %s', '%', i_model, model.name, model.ROIs)
    %     disp(model.lme)
    
    % results from some specific contrasts
    fprintf('%s %i - %s - %s\n', '%', i_model, model.name, model.ROIs)
    % effect of either linear or constant in either ROIs
    c = [...
        1 0 0 0 ;...
        0 1 0 0 ;...
        0 0 1 0 ;...
        0 0 0 1 ];
    message = 'effect of either linear or constant in either ROI';
    PVAL = test_and_print(model, c, pattern, message);
    
    if PVAL>.05
        c = [...
            1 0 0 0 ;...
            0 1 0 0];
        message = 'effect of either CST/LIN in ROI 1';
        PVAL = test_and_print(model, c, pattern, message);
        
        c = [...
            0 0 1 0 ;...
            0 0 0 1];
        message = 'effect of either CST/LIN in ROI 2';
        PVAL = test_and_print(model, c, pattern, message);
    else
        % effect of either CST/LIN
        c = [...
            1 0 1 0 ;...
            0 1 0 1 ];
        message = 'effect of either CST/LIN';
        PVAL = test_and_print(model, c, pattern, message);
        
        if PVAL>.05
            % effect of CST
            c = [1 0 1 0];
            message = 'effect of CST';
            PVAL = test_and_print(model, c, pattern, message);
            
            % effect of LIN
            c = [0 1 0 1];
            message = 'effect of LIN';
            PVAL = test_and_print(model, c, pattern, message);
        end
    end
    
    fprintf('\n\n\n')
end

end


function [PVAL,F,DF1,DF2] = test_and_print(model, c, pattern, message)
[PVAL,F,DF1,DF2] = coefTest(model.lme, c);
fprintf(pattern, ...
    message, ...
    DF1, DF2, ...
    F, PVAL);
end