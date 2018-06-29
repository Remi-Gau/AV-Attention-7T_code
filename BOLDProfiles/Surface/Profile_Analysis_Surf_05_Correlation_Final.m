%%
StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

addpath(genpath('/media/rxg243/BackUp2/AV_Integration_7T_2/SubFun'))

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

load(fullfile(StartFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'))

Subj2Include = true(13,1);
Subj2Include([4 11]) = false;

Xpos =  [1 1.5];

Legends1 = {'', 'Cst', '','','', '', '', '', '', 'Lin'};
Legends2 = {'ROI', 'mean', '(','SEM',')', 't value','p value', 'effect size', '', 'mean', '(','SEM',')', 't value','p value', 'effect size'};

COLOR_Subject = ColorSubject();

for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];

%%
close all

SubPlot = [1 3 2 4];

for iToPlot=2:3
    
    figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef\n(Fishcer trans)';
    end
    
    for iROI=1:1:numel(ROI)
        
        subplot(2,2,SubPlot(iROI))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCdt(iToPlot,1,iROI,:,Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        
        if isempty(ToPermute)
            [~,P] = ttest(Values2Plot);
        else
            Perms = ToPermute.*repmat(Values2Plot,[size(ToPermute,1),1]);
            Perms = mean(Perms,2);
            P = sum( abs( Perms-mean(Perms) ) > abs( mean(Values2Plot)-mean(Perms) ) ) / numel(Perms)
        end
        [H,P] = ttest(tmp');
        
        
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.5,sprintf(Sig));
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
        
%         h = plotSpread(tmp', ...
%             'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
%             'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
%         set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.5 0.5])
        
        t=title(ROI(iROI).name);
        set(t,'fontsize',12)
        
        t=ylabel(sprintf(YLabel));
        set(t,'fontsize',8)
    end
    
    mtit(['A vs V - ' suffix], 'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, fullfile(FigureFolder,'Baseline', ['Grp_Baseline_AvsV_' suffix '.tif']), '-dtiff')
    
end


SavedTxt = fullfile(FigureFolder,'Baseline','A_vs_V.csv');
fid = fopen (SavedTxt, 'w');
for iToPlot=2:3
    
    switch iToPlot
        case 2
            fprintf (fid, 'Slope');
        case 3
            fprintf (fid, 'CorrCoef');
    end
    
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    for iROI=1:numel(ROI)
        
        fprintf (fid, '%s,', ROI(iROI).name);
        
        tmp = squeeze(GrpBetaCdt(iToPlot,1,iROI,:,Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        
        for i=1:size(tmp,1)
            fprintf (fid, '%f,',nanmean(tmp(i,:)));
            fprintf (fid, '(,');
            fprintf (fid, '%f,',nansem(tmp(i,:)));
            fprintf (fid, '),');
            [~,P,~,STATS] = ttest(tmp(i,:), 0, 'alpha', 0.05);
            fprintf (fid, '%f,',STATS.tstat);
            fprintf (fid, '%f,',P);
            fprintf (fid, '%f,',abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
            fprintf (fid, ',');
        end
        
        fprintf (fid, '\n');
        
    end
    
    fprintf (fid, '\n\n\n');
end
fclose (fid);

%%
close all

SubPlot = [1 3 2 4];

Cdtion = [2 2 1 1];

CdtionName={'A vs AV-V','V vs AV-A'};

for iToPlot=2:3
    
    figure('name', 'CrossMod', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef\n(Fishcer trans)';
    end
    
    
    for iROI=1:1:numel(ROI)
        
        if iROI==1
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        elseif iROI==3
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        else
            SubPlotTitle = [];
        end
        
        subplot(2,2,SubPlot(iROI))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCrossSens(iToPlot,Cdtion(iROI),iROI,:,Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.5,sprintf(Sig));
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
        
%         h = plotSpread(tmp', ...
%             'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
%             'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
%         set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.25 .75])
        
        t=title(sprintf('%s\n%s',SubPlotTitle,ROI(iROI).name));
        set(t,'fontsize',12)
        
        t=ylabel(sprintf(YLabel));
        set(t,'fontsize',8)
        
        
    end
    
    mtit(suffix, 'fontsize', 14, 'xoff',0,'yoff',.05)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ['Grp_CrossSensory_' suffix '.tif']), '-dtiff')
end

SavedTxt = fullfile(FigureFolder,'CrossSensory','CrossSensory.csv');
fid = fopen (SavedTxt, 'w');
for iToPlot=2:3
    
    switch iToPlot
        case 2
            fprintf (fid, 'Slope');
        case 3
            fprintf (fid, 'CorrCoef');
    end
    
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    for iROI=1:numel(ROI)
        
        if iROI==1 || iROI==3
            fprintf (fid, '%s\n', CdtionName{Cdtion(iROI)});
        end
        
        fprintf (fid, '%s,', ROI(iROI).name);
        
        tmp = squeeze(GrpBetaCrossSens(iToPlot,Cdtion(iROI),iROI,:,Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        
        for i=1:size(tmp,1)
            fprintf (fid, '%f,',nanmean(tmp(i,:)));
            fprintf (fid, '(,');
            fprintf (fid, '%f,',nansem(tmp(i,:)));
            fprintf (fid, '),');
            [~,P,~,STATS] = ttest(tmp(i,:), 0, 'alpha', 0.05);
            fprintf (fid, '%f,',STATS.tstat);
            fprintf (fid, '%f,',P);
            fprintf (fid, '%f,',abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
            fprintf (fid, ',');
        end
        fprintf (fid, '\n');
        
    end
    
    fprintf (fid, '\n\n\n');
end


fclose (fid);



%% Cross vs Att
Name = {'AV-A'; 'AV-V'};
close all
CdtionName={'AV-A vs Att Mod (A-V)','AV-V vs Att Mod (V-A)'};

SubPlot = [1 3 2 4];

Cdtion = [1 1 2 2];

for iToPlot=2:3
    
    figure('name', 'CrossMod VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
            Ymin = -.25;
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef\n(Fishcer trans)';
            Ymin = -.25;
    end
    
    
    for iROI=1:1:numel(ROI)
        
        if iROI==1
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        elseif iROI==3
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        else
            SubPlotTitle = [];
        end
        
        subplot(2,2,SubPlot(iROI))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpCrossModAtt(iToPlot,Cdtion(iROI),iROI,:,Subj2Include));
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
            t = text(Xpos(i),0.6,sprintf(Sig));
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
        
%         h = plotSpread(tmp', ...
%             'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
%             'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
%         set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 Ymin 0.75])
        
        t=title(sprintf('%s\n%s',SubPlotTitle,ROI(iROI).name));
        set(t,'fontsize',12)
        
        t=ylabel(sprintf(YLabel));
        set(t,'fontsize',8)
        
    end
    
    mtit(suffix, 'fontsize', 14, 'xoff',0,'yoff',.05)
    
    
    print(gcf, fullfile(FigureFolder, 'CrossSensory', ['Grp_CrossMod_VS_Attention_' suffix '.tif']), '-dtiff')
    
end


SavedTxt = fullfile(FigureFolder,'CrossSensory','CrossSensory_VS_Att.csv');
fid = fopen (SavedTxt, 'w');
for iToPlot=2:3
    
    switch iToPlot
        case 2
            fprintf (fid, 'Slope');
        case 3
            fprintf (fid, 'CorrCoef');
    end
    
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    for iROI=1:numel(ROI)
        
        if iROI==1 || iROI==3
            fprintf (fid, '%s\n', CdtionName{Cdtion(iROI)});
        end
        
        fprintf (fid, '%s,', ROI(iROI).name);
        
        tmp = squeeze(GrpCrossModAtt(iToPlot,Cdtion(iROI),iROI,:,Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        if iROI>2
            tmp = tmp*-1;
        end
        
        for i=1:size(tmp,1)
            fprintf (fid, '%f,',nanmean(tmp(i,:)));
            fprintf (fid, '(,');
            fprintf (fid, '%f,',nansem(tmp(i,:)));
            fprintf (fid, '),');
            [~,P,~,STATS] = ttest(tmp(i,:), 0, 'alpha', 0.05);
            fprintf (fid, '%f,',STATS.tstat);
            fprintf (fid, '%f,',P);
            fprintf (fid, '%f,',abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
            fprintf (fid, ',');
        end
        fprintf (fid, '\n');
        
    end
    
    fprintf (fid, '\n\n\n');
end


fclose (fid);




%% Stim vs Att
close all

SubPlot = [1 3 2 4];

CdtionName={'A stim vs Att Mod','V stim vs Att Mod','AV Stim vs AttMod'};

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
            
%             h = plotSpread(tmp', ...
%                 'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
%                 'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
%             set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
            
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


SavedTxt = fullfile(FigureFolder,'Attention','Stim_VS_Attention.csv');
fid = fopen (SavedTxt, 'w');

for iCdt=1:numel(CdtionName)
    
    fprintf (fid, '%s\n', CdtionName{iCdt});
    
    for iToPlot=2:3
        
        switch iToPlot
            case 2
                fprintf (fid, 'Slope');
            case 3
                fprintf (fid, 'CorrCoef');
        end
        
        for i=1:length(Legends1)
            fprintf (fid, '%s,', Legends1{i});
        end
        fprintf (fid, '\n');
        for i=1:length(Legends2)
            fprintf (fid, '%s,', Legends2{i});
        end
        fprintf (fid, '\n');
        
        for iROI=1:numel(ROI)
            
            fprintf (fid, '%s,', ROI(iROI).name);
            
            tmp = squeeze(GrpBetaStimAtt(iToPlot,iCdt,iROI,:,Subj2Include));
            if iToPlot==3
                tmp = atanh(tmp);
            end
            if iROI>2
                tmp = tmp*-1;
            end
            
            for i=1:size(tmp,1)
                fprintf (fid, '%f,',nanmean(tmp(i,:)));
                fprintf (fid, '(,');
                fprintf (fid, '%f,',nansem(tmp(i,:)));
                fprintf (fid, '),');
                [~,P,~,STATS] = ttest(tmp(i,:), 0, 'alpha', 0.05);
                fprintf (fid, '%f,',STATS.tstat);
                fprintf (fid, '%f,',P);
                fprintf (fid, '%f,',abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
                fprintf (fid, ',');
            end
            fprintf (fid, '\n');
            
        end
        
        fprintf (fid, '\n\n');
        
    end
    
    fprintf (fid, '\n\n\n');
end


fclose (fid);


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
            
%             h = plotSpread(tmp', ...
%                 'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
%                 'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
%             set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
            
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




SavedTxt = fullfile(FigureFolder,'Attention','StimSpe_VS_Attention.csv');
fid = fopen (SavedTxt, 'w');

for iCdt=1:numel(CdtionName)
    
    fprintf (fid, '%s\n', CdtionName{iCdt});
    
    for iToPlot=2:3
        
        switch iToPlot
            case 2
                fprintf (fid, 'Slope');
            case 3
                fprintf (fid, 'CorrCoef');
        end
        
        for i=1:length(Legends1)
            fprintf (fid, '%s,', Legends1{i});
        end
        fprintf (fid, '\n');
        for i=1:length(Legends2)
            fprintf (fid, '%s,', Legends2{i});
        end
        fprintf (fid, '\n');
        
        for iROI=1:numel(ROI)
            
            fprintf (fid, '%s,', ROI(iROI).name);
            
            tmp = squeeze(GrpBetaStimSpeAtt(iToPlot,iCdt,iROI,:,Subj2Include));
            if iToPlot==3
                tmp = atanh(tmp);
            end
            if iROI>2
                tmp = tmp*-1;
            end
            
            for i=1:size(tmp,1)
                fprintf (fid, '%f,',nanmean(tmp(i,:)));
                fprintf (fid, '(,');
                fprintf (fid, '%f,',nansem(tmp(i,:)));
                fprintf (fid, '),');
                [~,P,~,STATS] = ttest(tmp(i,:), 0, 'alpha', 0.05);
                fprintf (fid, '%f,',STATS.tstat);
                fprintf (fid, '%f,',P);
                fprintf (fid, '%f,',abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
                fprintf (fid, ',');
            end
            fprintf (fid, '\n');
            
        end
        
        fprintf (fid, '\n\n');
        
    end
    
    fprintf (fid, '\n\n\n');
end


fclose (fid);
