%%
clc; clear;

CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))

FigureFolder = fullfile(CodeFolder,'Figures');

SourceFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

load(fullfile(SourceFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'))

Subj2Include = true(11,1);

Xpos =  [1 1.5];


%%
close all
CdtionName={'A vs V','A vs AV','V vs AV'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCdt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        for i=1:size(tmp,1)
            [P,H] = SignPermTest(tmp(i,:));
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),1.75,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H==1
                set(t,'color','r');
            end
        end
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.5 1.75])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',12)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

% print(gcf, fullfile(FigureFolder,'Baseline', 'Grp_CorrBaseline.tif'), '-dtiff')

%%
close all
CdtionName={'A vs AV-V','V vs AV-A'};
figure('name', 'CrossSensory', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCrossSens(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),1.5,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.25 1.5])
    end
end

subplot(numel(ROI),2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),2,3)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,5)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,7)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder,'CrossSensory', 'Grp__CorrCrossSensory.tif'), '-dtiff')


%%
close all
CdtionName={'Att Mod_{A stim} vs Att Mod_{V stim}','Att Mod_{A stim} vs Att Mod_{AV stim}','Att Mod_{V stim} vs AttMod_{AV stim}'};
figure('name', 'Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.55,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.20 0.55])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_CorrAttention.tif'), '-dtiff')


%% Stim vs Att
close all
CdtionName={'A stim vs Att Mod','V stim vs Att Mod','AV Stim vs AttMod'};
figure('name', 'Stiim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaStimAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_Corr_Stim_VS_Attention.tif'), '-dtiff')


%% Stim vs stim specific Att
close all
CdtionName={'A stim vs Att Mod A stim','V stim vs Att Mod V stim','AV Stim vs AttMod AV stim'};
figure('name', 'Stim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaStimSpeAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_Corr_StimSpe_VS_Attention.tif'), '-dtiff')



%% Cross vs Att
Name = {'AV-A'; 'AV-V'};
close all
CdtionName={'AV-A vs Att Mod','AV-V vs Att Mod'};
figure('name', 'CrossMod VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpCrossModAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),2,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),2,3)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,5)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'CrossSensory', 'Grp_Corr_CrossMod_VS_Attention.tif'), '-dtiff')
