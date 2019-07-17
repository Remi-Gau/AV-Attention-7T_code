function [model, Y_legend] = linear_mixed_model(data_rois)
% Runs a F-test on for the shape parameters (encoded in 2 different regressors)

nb_rois = numel(data_rois);

%% extract data
Y = [];
X = [];

Y_legend = {};

for iROI = 1:nb_rois
    
    DATA = data_rois{iROI};
    
    nb_subjects = size(DATA.Betas,1);
    
    Y_legend = gen_y_legend(Y_legend, DATA, nb_subjects);
    
    % concatenate data over ROIs
    y = DATA.Betas;
    y = y(:);
    Y = [Y ; y]; %#ok<*AGROW>
    
    % design matrix:
    % reg 1 : cst / lin
    % reg 2 : ROI
    % reg 3 : subject
    X = [ X; ...
        ones(nb_subjects,1) ones(nb_subjects,1)*iROI [1:nb_subjects]' ; ...
        ones(nb_subjects,1)*2 ones(nb_subjects,1)*iROI [1:nb_subjects]' ];
    
end

%% plot
% figure(1)
% colormap gray
%
% subplot(1,2,1)
% imagesc(Y)
% set(gca, 'ytick', 1:numel(Y), 'yticklabel', Y_legend)
% title('data')
%
% subplot(1,2,2)
% imagesc(X)
% set(gca, 'xtick', 1:4, 'xticklabel', {'cst/lin'; 'ROI'; 'subject'})
% title('design matrix')


%% LMM
% subject is random
% lin/const and ROI are fixed

Subject = X(:,3); % Fixed effect predictor: 'ShapeParameter','ROI','Intercept'
G = Subject; % RandomEffect Group
Z = ones(size(X,1),1); % random effect predictor

x = [X(:,1:2) ones(size(X,1),1)];

lme = fitlmematrix(x, Y, Z, G,...
    'FixedEffectPredictors',...
        {'ShapeParameter','ROI','Intercept'},...
    'RandomEffectPredictors',...
        {{'Intercept'}},...
    'RandomEffectGroups',...
        {'Subject'});

model.X = x;
model.G = G;
model.Z = Z;
model.Y = Y;
model.lme = lme;

end

