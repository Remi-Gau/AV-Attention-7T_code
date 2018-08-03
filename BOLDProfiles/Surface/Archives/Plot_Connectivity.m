%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

COLOR_Subject = ColorSubject();

NbLayers = 7;
NbLayers = NbLayers+1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigureFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Figures\Profiles\Surfaces\Baseline';

Results_Folder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Results\Profiles\Surfaces';
load(fullfile(Results_Folder,strcat('Connectivity_Surf_', num2str(NbLayers), '_layers.mat')))


Colors = [...
    1 0 0;
    0 1 0;
    0 0 1;
    .5 0 0;
    0 .5 0;
    0 0 .5];
FigDim = [100 100 1000 700];

close all

%% run permutation test on A1 cst
sets = {};
for iSub=1:11
    sets{iSub} = [-1 1]; %#ok<*AGROW>
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


% %% Plot slope for just V-Fix averaged over attention and based on average value at each vertex
% figure('name', 'rank', 'Position', FigDim, 'Color', [1 1 1]);
%
% hold on
% grid on
%
% subj2incl = logical([0 0 1 1 1 1 1 1 0 1 1]);
%
% for iSub = 1:size(B_rank,2)
%     plot(B_rank_avg(2:end,iSub,2),V1_Act_rank(iSub,2:end,2),'color',...
%         COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
% end
%
% plot(mean(B_rank_avg(2:end,subj2incl,2),2),mean(V1_Act_rank(:,2:end,2)),'linewidth',1.5, ...
%     'color', 'b')


%% Plot slope for each of the 6 conditions and based on Cst S param at each vertex
% this is always plotted against the value for V-Fix averaged over
% attention (based on Cst S param at each vertex)
close all
for Eig_or_Mean = 1:2
    
    figure('name', 'rank', 'Position', FigDim, 'Color', [1 1 1]);
    
    if Eig_or_Mean==1
        Suffix='Eig';
    else
        Suffix='Mean' ;
    end
    
    for iCdt = 1:6
        
        subplot(3,3,iCdt)
        
        hold on
        grid on
        
        for iSub = 1:size(B_rank,2)
            plot(mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),mean(V1_Act_rank(iSub,2:end,2),4),'color',...
                COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(V1_Act_rank(:,2:end)),....
        %             nansem(V1_Act_rank(:,2:end))*0,nansem(V1_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),mean(mean(V1_Act_rank(:,2:end,2),4)),'linewidth',1.5, ...
            'color', Colors(iCdt,:))
        
        t=herrorbar(mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),16,nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),14,nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        title('Connectivity ; slope')

            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 14 16],  'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'})
            axis([-.05 .05 -4 16.5])

        
        if iCdt==1 || iCdt==4
            ylabel('[V-Fix]')
        end
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    subplot(3,3,7:9)
    hold on
    grid on
    t=plot(repmat(mean(V1_Act_rank(:,2:end,2))',1,6),squeeze(mean(B_rank(2:end,:,:,2,Eig_or_Mean),2)));
    for i=1:numel(t)
        set(t(i), 'linewidth', 1, 'markersize', 1, 'color', Colors(i,:))
    end
    plot([-50 50],[0 0], '--k')
    plot([0 0],[-50 50], '--k')
    
    xlabel('[V-Fix]')
    ylabel('Connectivity ; slope')
    
    axis([-4 12 -.01 0.01])

    mtit(['Connectvity V1 to A1 ' Suffix],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectvity_V1_to_A1_CV_' Suffix ' .tif']), '-dtiff')
    
end


%% Plot slope for each of the 2 X 2 conditions and their slope average
% this is always plotted against the value for V-Fix averaged over
% attention (based on Cst S param at each vertex)
for Eig_or_Mean = 1:2
    
    figure('name', 'rank', 'Position', FigDim, 'Color', [1 1 1]);
    
    if Eig_or_Mean==1
        Suffix='Eig';
    else
        Suffix='Mean' ;
    end
    
    for iCdt = [1 2 4 5]
        
        subplot(3,3,iCdt)
        
        hold on
        grid on
        
        for iSub = 1:size(B_rank,2)
            plot(mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),mean(V1_Act_rank(iSub,2:end,2),4),'color',...
                COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(V1_Act_rank(:,2:end)),....
        %             nansem(V1_Act_rank(:,2:end))*0,nansem(V1_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),mean(mean(V1_Act_rank(:,2:end,2),4)),'linewidth',1.5, ...
            'color', Colors(iCdt,:))
        
        t=herrorbar(mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),16,nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),14,nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
            'ytick', [-8:4:12 14 16],  'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'})
        axis([-.05 .05 -4 16.5])

        if iCdt==1 || iCdt==4
            ylabel('[V-Fix]')
        end
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    for iCdt = 1:4
        
        switch iCdt
            case 1
                subplot(3,3,3)
                CdtToAvg = [1 2];
            case 2
                subplot(3,3,6)
                CdtToAvg = [4 5];
            case 3
                subplot(3,3,7)
                CdtToAvg = [1 4];
            case 4
                subplot(3,3,8)
                CdtToAvg = [2 5];
        end
        
        hold on
        grid on
        
        for iSub = 1:size(B_rank,2)
            plot(mean(mean(B_rank(2:end,iSub,CdtToAvg,2,Eig_or_Mean),6),3),...
                mean(V1_Act_rank(iSub,2:end,2),4),'color',...
                COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        end
        
        %         t=errorbar(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),mean(V1_Act_rank(:,2:end)),....
        %             nansem(V1_Act_rank(:,2:end))*0,nansem(V1_Act_rank(:,2:end))*0,...
        %             nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1.5, 'markersize', 2,'color', 'k')
        
        plot(mean(mean(mean(B_rank(2:end,:,CdtToAvg,2,Eig_or_Mean),6),3),2),...
            mean(mean(V1_Act_rank(:,2:end,2),4)),'linewidth',1.5, ...
            'color', 'k')
        
        t=herrorbar(mean(mean(mean(B_act(:,CdtToAvg,2,Eig_or_Mean),5),2)),16,...
            nansem(mean(mean(B_act(:,CdtToAvg,2,Eig_or_Mean),5),2)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', 'k')
        t=herrorbar(mean(mean(mean(B_deact(:,CdtToAvg,2,Eig_or_Mean),5),2)),14,...
            nansem(mean(mean(B_deact(:,CdtToAvg,2,Eig_or_Mean),5),2)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', 'k')
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        % compute perm test
        tmp = [mean(mean(B_deact(:,CdtToAvg,2,Eig_or_Mean),5),2) ...
            mean(mean(B_act(:,CdtToAvg,2,Eig_or_Mean),5),2)];
        if iCdt==1 || iCdt==3
            tmp2 = diff(tmp,[],2);
        else
            tmp2 = [tmp2 diff(tmp,[],2)];
            
            tmp2 = diff(tmp2,[],2);
            
            for iPerm = 1:size(ToPermute,1)
                tmp3 = ToPermute(iPerm,:);
                tmp3 = repmat(tmp3',1,size(tmp2,2));
                Perms(iPerm,:) = mean(tmp2.*tmp3);  %#ok<*AGROW>
            end
            P = sum( ...
                abs( Perms ) > ...
                repmat( abs(mean(tmp2)), size(Perms,1),1)  ) ...
                / size(Perms,1);
            
            Sig = [];
            if P<0.001
                Sig = 'p<0.001';
            else
                Sig = sprintf('p=%.3f',P);
            end
            t = text(.01,16,sprintf(Sig));
            set(t,'fontsize',10);
            
            if P<0.05
                set(t,'fontsize',Fontsize,'fontweight','bold')
                set(t,'color','r');
            end
        end

        set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
            'ytick', [-8:4:12 14 16],  'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'})
        axis([-.05 .05 -4 16.5])
        
        switch iCdt
            case 1
                title(sprintf('Connectivity ; slope\nA attention'))
            case 2
                title(sprintf('Connectivity ; slope\nV attention'))
            case 3
                title(sprintf('Connectivity ; slope\nA stim'))
            case 4
                title(sprintf('Connectivity ; slope\nV stim'))
        end
        
    end
    
    mtit(['Connectvity V1 to A1 ' Suffix ' - slope averaging'],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectvity_V1_to_A1_CV_slope_averaging' Suffix ' .tif']), '-dtiff')
    
end

%% Plot slope for each of the 2 X 2 conditions and the slope resulting of their concatenation
% this is always plotted against the value for V-Fix averaged over
% attention (based on Cst S param at each vertex)
for Eig_or_Mean = 1:2
    
    figure('name', 'rank', 'Position', FigDim, 'Color', [1 1 1]);
    
    if Eig_or_Mean==1
        Suffix='Eig';
    else
        Suffix='Mean' ;
    end
    
    for iCdt = [1 2 4 5]
        
        subplot(3,3,iCdt)
        
        hold on
        grid on
        
        for iSub = 1:size(B_rank,2)
            plot(mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),mean(V1_Act_rank(iSub,2:end,2),4),'color',...
                COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(V1_Act_rank(:,2:end)),....
        %             nansem(V1_Act_rank(:,2:end))*0,nansem(V1_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),mean(mean(V1_Act_rank(:,2:end,2),4)),'linewidth',1.5, ...
            'color', Colors(iCdt,:))
        
        t=herrorbar(mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),16,nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),14,nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')

        set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
            'ytick', [-8:4:12 14 16],  'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'})
        axis([-.05 .05 -4 16.5])
        
        if iCdt==1 || iCdt==4
            ylabel('[V-Fix]')
        end
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    for iCdt = 1:4
        
        switch iCdt
            case 1
                subplot(3,3,3)
                RankData = B_rank_att;
                ActData = B_act_att;
                DeactData = B_deact_att;
                CdtToPlot = 1;
            case 2
                subplot(3,3,6)
                RankData = B_rank_att;
                ActData = B_act_att;
                DeactData = B_deact_att;
                CdtToPlot = 2;
            case 3
                subplot(3,3,7)
                RankData = B_rank_sens;
                ActData = B_act_sens;
                DeactData = B_deact_sens;
                CdtToPlot = 1;
            case 4
                subplot(3,3,8)
                RankData = B_rank_sens;
                ActData = B_act_sens;
                DeactData = B_deact_sens;
                CdtToPlot = 2;
        end
        
        hold on
        grid on
        
        for iSub = 1:size(RankData,2)
            plot(mean(RankData(2:end,iSub,CdtToPlot,2,Eig_or_Mean),6),mean(V1_Act_rank(iSub,2:end,2),4),'color',...
                COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        end
        
        %         t=errorbar(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),mean(V1_Act_rank(:,2:end)),....
        %             nansem(V1_Act_rank(:,2:end))*0,nansem(V1_Act_rank(:,2:end))*0,...
        %             nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1.5, 'markersize', 2,'color', 'k')
        
        plot(mean(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),6),2),mean(mean(V1_Act_rank(:,2:end,2),4)),'linewidth',1.5, ...
            'color', 'k')
        
        t=herrorbar(mean(mean(ActData(:,CdtToPlot,2,Eig_or_Mean),5)),16,nansem(mean(ActData(:,CdtToPlot,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', 'k')
        t=herrorbar(mean(mean(DeactData(:,CdtToPlot,2,Eig_or_Mean),5)),14,nansem(mean(DeactData(:,CdtToPlot,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', 'k')
        plot([-15 15],[0 0], '--k')
        plot([0 0],[-15 15], '--k')
        
        
        % compute perm test
        tmp = [mean(DeactData(:,CdtToPlot,2,Eig_or_Mean),5) ...
            mean(ActData(:,CdtToPlot,2,Eig_or_Mean),5)];
        if iCdt==1 || iCdt==3
            tmp2 = diff(tmp,[],2);
        else
            tmp2 = [tmp2 diff(tmp,[],2)];
            
            tmp2 = diff(tmp2,[],2);
            
            for iPerm = 1:size(ToPermute,1)
                tmp3 = ToPermute(iPerm,:);
                tmp3 = repmat(tmp3',1,size(tmp2,2));
                Perms(iPerm,:) = mean(tmp2.*tmp3);  %#ok<*AGROW>
            end
            P = sum( ...
                abs( Perms ) > ...
                repmat( abs(mean(tmp2)), size(Perms,1),1)  ) ...
                / size(Perms,1);
            
            Sig = [];
            if P<0.001
                Sig = 'p<0.001';
            else
                Sig = sprintf('p=%.3f',P);
            end
            t = text(.01,16,sprintf(Sig));
            set(t,'fontsize',10);
            
            if P<0.05
                set(t,'fontsize',10,'fontweight','bold')
                set(t,'color','r');
            end
        end


        set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
            'ytick', [-8:4:12 14 16],  'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'})
        axis([-.08 .08 -4 16.5])

        switch iCdt
            case 1
                title(sprintf('Connectivity ; slope\nA attention'))
            case 2
                title(sprintf('Connectivity ; slope\nV attention'))
            case 3
                title(sprintf('Connectivity ; slope\nA stim'))
            case 4
                title(sprintf('Connectivity ; slope\nV stim'))
        end
        
    end
    
    mtit(['Connectvity V1 to A1 ' Suffix ' - concatenation'],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectvity_V1_to_A1_CV_concat' Suffix ' .tif']), '-dtiff')
    
end
