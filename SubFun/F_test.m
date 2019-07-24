function [F_semi_partial_coef, p_semi_partial_coef, df, dferror] = F_test(data_rois)
% Runs a F-test on for the shape parameters (encoded in 2 different regressors)

% Code massively inspired from Cyril's website
% http://www.sbirc.ed.ac.uk/cyril/glm/GLM_lectures.html

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
    % reg 1 : cst
    % reg 2 : lin
    % reg 3:end : subject specific regressors
    X = [ X; ...
        ones(nb_subjects,1) zeros(nb_subjects,1) eye(nb_subjects); ...
        zeros(nb_subjects,1) ones(nb_subjects,1) eye(nb_subjects)];
    
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
% set(gca, 'xtick', [1, 2], 'xticklabel', {'cst'; 'lin'})
% title('design matrix')



%% do GLM and stats
[B, Yhat, Res, SStotal, SSeffect, SSerror, R2, df, dferror, F, p] = do_glm(X, Y); %#ok<*ASGLU>

Xreduced    = X(:,3:end); % reduced model all minus 2 1st regressor
[B_red, Yhat_red, Res_red, SStotal_red, SSeffect_red, SSerror_red, R2_red, df_red, dferror_red, F_red, p_red] =...
    do_glm(Xreduced, Y);

% Compare reduced model and full model to see how much variance is
% explained by Cst and Lin regressors
Semi_Partial_corr_coef = R2 - R2_red;
dfe_semi_partial_coef  = df - df_red;
F_semi_partial_coef    = (Semi_Partial_corr_coef*dferror) / ...  % variance explained by x1
    ((1-R2)*dfe_semi_partial_coef); % unexplained variance overall
p_semi_partial_coef    = 1 - fcdf(Semi_Partial_corr_coef, df, dfe_semi_partial_coef); % note df is from the full model

end

function [B, Yhat, Res, SStotal, SSeffect, SSerror, R2, df, dferror, F, p] = do_glm(X, Y)
% Runs GLM and gets F and p value for the whole model.

B    = inv(X'*X)*X'*Y;
Yhat = X*B;
Res  = Y - Yhat;

SStotal = norm(Y - mean(Y)).^2;
SSeffect = norm(Yhat - mean(Yhat)).^2;
SSerror  = norm(Res-mean(Res)).^2;
R2      = SSeffect / SStotal;

df      = rank(X)-1;
dferror = length(Y) - df - 1;
F       = (SSeffect / df) / (SSerror / dferror);
p       = 1 - fcdf(F,df,dferror);

end