function [F_semi_partial_coef, p_semi_partial_coef, df, dferror] = F_test(data_rois)

% http://www.sbirc.ed.ac.uk/cyril/glm/GLM_lectures.html

nb_rois = numel(data_rois);

%% extract data
Y = [];
X = [];

Y_legend = {};

for iROI = 1:nb_rois
    
    DATA = data_rois{iROI};
    
    nb_subjects = size(DATA.Betas,1);
    
    Y_legend{end+1} = [DATA.Name ' ; sub-01 ; CST' ]; %#ok<*AGROW>
    for i=1:nb_subjects-2
        Y_legend{end+1} = '...';
    end
    Y_legend{end+1} = [DATA.Name ' ; sub-11 ; CST' ];
    Y_legend{end+1} = [DATA.Name ' ; sub-01 ; LIN' ];
    for i=1:nb_subjects-2
        Y_legend{end+1} = '...';
    end
    Y_legend{end+1} = [DATA.Name ' ; sub-11 ; LIN' ];
    
    y = DATA.Betas;
    y = y(:);
    
    Y = [Y ; y];
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
[B, Yhat, Res, SStotal, SSeffect, SSerror, R2, df, dferror, F, p] = do_glm(X, Y);

Xreduced    = X(:,3:end); % reduced model all minus 2 1st regressor
[B_red, Yhat_red, Res_red, SStotal_red, SSeffect_red, SSerror_red, R2_red, df_red, dferror_red, F_red, p_red] =...
    do_glm(Xreduced, Y);

Semi_Partial_corr_coef = R2 - R2_red;
dfe_semi_partial_coef  = df - df_red;
F_semi_partial_coef    = (Semi_Partial_corr_coef*dferror) / ...  % variance explained by x1
    ((1-R2)*dfe_semi_partial_coef); % unexplained variance overall
p_semi_partial_coef    = 1 - fcdf(Semi_Partial_corr_coef, df, dfe_semi_partial_coef); % note df is from the full model

end

function [B, Yhat, Res, SStotal, SSeffect, SSerror, R2, df, dferror, F, p] = do_glm(X, Y)

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