% plots some specific vertex wise correlation and also the slope of the
% linear relationship between the 2 variables

% corellation coefficient: fisher transformation applied
% significance computed with exact sign permutation test

% More complex figures still need further refactoring

%%
clc; clear; close all

CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))
Get_dependencies('/home/remi')

SourceFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

load(fullfile(SourceFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'))

opt.FigDim = [100, 100, 1000, 1500];
opt.Visibility = 'on';
opt.Xpos =  [1 1.5];
opt.Subj2Include = true(11,1);
opt.fontsize = 10;
opt.FigureFolder = fullfile(CodeFolder,'Figures');
opt.print = 0;


%%
opt.name = 'Grp_Baseline_AvsV';
Cdtion = [1 1 1 1];
CdtionName= {'', ''};
plot_group_vertex_slope_corr(GrpBetaCdt, CdtionName, Cdtion, ROI, opt)


opt.name = 'Grp_CrossSensory';
Cdtion = [2 2 1 1];
CdtionName={'A vs AV-V','V vs AV-A'};
plot_group_vertex_slope_corr(GrpBetaCrossSens, CdtionName, Cdtion, ROI, opt)


%% Cross vs Att
opt.name = 'Grp_CrossMod_VS_Attention';
CdtionName={'AV-A vs Att Mod (A-V)','AV-V vs Att Mod (V-A)'};
Cdtion = [1 1 2 2];

plot_group_vertex_slope_corr(GrpCrossModAtt, CdtionName, Cdtion, ROI, opt)





%% Stim vs Att



CdtionName={'A stim vs Att Mod','V stim vs Att Mod','AV Stim vs AttMod'};

plot_group_vertex_slope_corr(GrpCrossModAtt, CdtionName, Cdtion, ROI, opt)

for iToPlot=2:3
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
            Ymin = -.35;
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef';
            Ymin = -.35;
    end
    
    figure('name', 'Stiim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
    
    for iCdt=1:numel(CdtionName)
        for iROI=1:1:numel(ROI)
            
            subplot(numel(ROI),3,iCdt+3*(iROI-1))
            hold on
            
            plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
            
            tmp = squeeze(GrpBetaStimAtt(iToPlot,iCdt,iROI,:,Subj2Include));
            if iToPlot==3
                tmp = atanh(tmp);
            end
            if iROI>2
                tmp = tmp*-1;
            end
            [H,P] = ttest(tmp');
            for i=1:numel(P)
                if P(i)<0.001
                    Sig = sprintf('\np<0.001 ');
                else
                    Sig = sprintf('\np=%.3f ',P(i));
                end
                t = text(Xpos(i),0.7,sprintf(Sig));
                set(t,'fontsize',8);
                
                if H(i)==1
                    set(t,'color','r');
                end
            end
            
            h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
            set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
            
            for i=1:size(tmp,2)
                plot(1.025+i*.025,tmp(1,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
                plot(1.525+i*.025,tmp(2,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
            end
            

            
            set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
                'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
                'ygrid', 'on')
            axis([0.9 1.8 Ymin 0.8])
        end
    end
    
    subplot(numel(ROI),3,1)
    t=title(CdtionName{1});
    set(t,'fontsize',8)
    t=ylabel(sprintf('%s \n %s',ROI(1).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,2)
    t=title(CdtionName{2});
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,3)
    t=title(CdtionName{3});
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,4)
    t=ylabel(sprintf('%s \n %s',ROI(2).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,7)
    t=ylabel(sprintf('%s \n %s',ROI(3).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,10)
    t=ylabel(sprintf('%s \n %s',ROI(4).name, YLabel));
    set(t,'fontsize',8)
    
    mtit(suffix, 'fontsize', 14, 'xoff',0,'yoff',.05)
    
    print(gcf, fullfile(FigureFolder, 'Attention', ['Grp_Stim_VS_Attention_' suffix '.tif']), '-dtiff')
    
end


%% Stim vs stim specific Att
close all

SubPlot = [1 3 2 4];

CdtionName={'A stim vs Att Mod A stim','V stim vs Att Mod V stim','AV Stim vs AttMod AV stim'};
for iToPlot=2:3
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
            Ymin = -.5;
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef';
            Ymin = -.75;
    end
    
    figure('name', 'Stim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
    
    for iCdt=1:numel(CdtionName)
        for iROI=1:1:numel(ROI)
            
            subplot(numel(ROI),3,iCdt+3*(iROI-1))
            hold on
            
            plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
            
            tmp = squeeze(GrpBetaStimSpeAtt(iToPlot,iCdt,iROI,:,Subj2Include));
            if iToPlot==3
                tmp = atanh(tmp);
            end
            if iROI>2
                tmp = tmp*-1;
            end
            [H,P] = ttest(tmp');
            for i=1:numel(P)
                if P(i)<0.001
                    Sig = sprintf('\np<0.001 ');
                else
                    Sig = sprintf('\np=%.3f ',P(i));
                end
                t = text(Xpos(i),0.75,sprintf(Sig));
                set(t,'fontsize',8);
                
                if H(i)==1
                    set(t,'color','r');
                end
            end
            
            h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
            set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
            
            for i=1:size(tmp,2)
                plot(1.025+i*.025,tmp(1,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
                plot(1.525+i*.025,tmp(2,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
            end

            
            set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
                'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
                'ygrid', 'on')
            axis([0.9 1.8 -0.35 0.75])
        end
    end
    
    subplot(numel(ROI),3,1)
    t=title(CdtionName{1});
    set(t,'fontsize',8)
    t=ylabel(sprintf('%s \n %s',ROI(1).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,2)
    t=title(CdtionName{2});
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,3)
    t=title(CdtionName{3});
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,4)
    t=ylabel(sprintf('%s \n %s',ROI(2).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,7)
    t=ylabel(sprintf('%s \n %s',ROI(3).name, YLabel));
    set(t,'fontsize',8)
    
    subplot(numel(ROI),3,10)
    t=ylabel(sprintf('%s \n %s',ROI(4).name, YLabel));
    set(t,'fontsize',8)
    
    mtit(suffix, 'fontsize', 14, 'xoff',0,'yoff',.05)
    
    print(gcf, fullfile(FigureFolder, 'Attention', ['Grp_StimSpe_VS_Attention_' suffix '.tif']), '-dtiff')
    
    
end