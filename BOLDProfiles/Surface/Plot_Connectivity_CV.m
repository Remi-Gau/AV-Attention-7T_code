%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

COLOR_Subject = ColorSubject();

V1_2_A1 = 0;
if V1_2_A1==1
    ROI_Suffix = 'V1_2_A1';
    Cdt_range_to_plot = 5;
    Y_all = 19;
    Y_act = 16;
    Y_deact = 14;
else
    ROI_Suffix = 'A1_2_V1';
    Cdt_range_to_plot = 4;
    Y_act = 10;
    Y_deact = 8;
    Y_all = 13;
end

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
load(fullfile(Results_Folder,strcat('Connectivity_', ROI_Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')))

% load(fullfile(Results_Folder,strcat('Connectivity_all_', ROI_Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')))

Colors = [...
    1 0 0;
    0 1 0;
    0 0 1;
    .5 0 0;
    0 .5 0;
    0 0 .5];
FigDim = [50 50 1200 750];

close all

%% run permutation test on A1 cst
sets = {};
for iSub=1:11
    sets{iSub} = [-1 1]; %#ok<*AGROW>
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


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
        
        %         for iSub = 1:size(B_rank,2)
        %             plot(...
        %                 mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),...
        %                 mean(Seed_Act_rank(iSub,2:end,Cdt_range_to_plot),4),...
        %                 'color',COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        %         end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(Seed_Act_rank(:,2:end)),....
        %             nansem(Seed_Act_rank(:,2:end))*0,nansem(Seed_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(...
            mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),...
            mean(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot),4)),...
            'linewidth',1.5, 'color', Colors(iCdt,:))
        
        t=herrorbar(...
            mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),...
            Y_act,...
            nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),...
            Y_deact,...
            nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        title('Connectivity ; slope')
        
        if V1_2_A1==1
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -4 19.5])
        else
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-2:2:6 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-2';'0';'2';'4';'6'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -2 13.5])
        end
        
        if iCdt==1 || iCdt==4
            if V1_2_A1
                ylabel('[V-Fix]')
            else
                ylabel('[A-Fix]')
            end
        end
        
        % compute perm test Y_act VS Y deact
        tmp = diff([mean(B_act(:,iCdt,2,Eig_or_Mean),5) ...
            mean(B_deact(:,iCdt,2,Eig_or_Mean),5)],[],2);
        
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = 'p_{act / deact}<0.001';
        else
            Sig = sprintf('p_{act / deact}=%.3f',P);
        end
        t = text(.025,Y_act,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp2 tmp3
        
        
        % compute perm test Y_all
        tmp = mean(B_all(:,iCdt,2,Eig_or_Mean),5);
        
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = 'p_{all}<0.001';
        else
            Sig = sprintf('p_{all}=%.3f',P);
        end
        t = text(.025,Y_all,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp2 tmp3
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    subplot(3,3,7:9)
    hold on
    grid on
    t=plot(...
        repmat(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot))',1,6),...
        squeeze(mean(B_rank(2:end,:,:,2,Eig_or_Mean),2)));
    for i=1:numel(t)
        set(t(i), 'linewidth', 1, 'markersize', 1, 'color', Colors(i,:))
    end
    plot([-50 50],[0 0], '--k')
    plot([0 0],[-50 50], '--k')
    
    xlabel('[A-Fix]')
    ylabel('Connectivity ; slope')
    if V1_2_A1==1
        axis([-2 4 0 0.035])
    else
        axis([-2 4 -0.005 0.02])
    end
    
    mtit(['Connectivity ' strrep(ROI_Suffix,'_2_',' to ') ' - ' Suffix],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectivity_' ROI_Suffix '_' Suffix ' .tif']), '-dtiff')
    
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
        
        %         for iSub = 1:size(B_rank,2)
        %             plot(...
        %                 mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),...
        %                 mean(Seed_Act_rank(iSub,2:end,Cdt_range_to_plot),4),...
        %                 'color',COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        %         end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(Seed_Act_rank(:,2:end)),....
        %             nansem(Seed_Act_rank(:,2:end))*0,nansem(Seed_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(...
            mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),...
            mean(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot),4)),...
            'linewidth',1.5,'color', Colors(iCdt,:))
        
        t=herrorbar(...
            mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),...
            Y_act,...
            nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),...
            Y_deact,...
            nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        List_B(:,iCdt) = mean(B_all(:,iCdt,2,Eig_or_Mean),5);
        
        if V1_2_A1==1
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -4 19.5])
        else
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-2:2:6 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-2';'0';'2';'4';'6'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -2 13.5])
        end
        
        if iCdt==1 || iCdt==4
            if V1_2_A1
                ylabel('[V-Fix]')
            else
                ylabel('[A-Fix]')
            end
        end
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    
    %     i. Att A > Att V
    %     ii.  Stim A > V
    %     iii. Stim A > V for AttA > Att V
    List_B(:,3) = [];
    
    for i=1:3
        switch i
            case 1
                subplot(3,3,2)
                P_suffix = 'Att_A > Att_V';
                tmp = sum(List_B(:,[1 2]),2)-sum(List_B(:,[3 4]),2);
            case 2
                subplot(3,3,4)
                P_suffix = 'Stim A > V';
                tmp = sum(List_B(:,[1 3]),2)-sum(List_B(:,[2 4]),2);
            case 3
                subplot(3,3,5)
                P_suffix = 'Stim A > V for Att_A > Att_V';
                tmp = diff(List_B(:,[1 3]),[],2)-diff(List_B(:,[2 4]),[],2);
        end
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = sprintf('p_{%s - all}<0.001',P_suffix);
        else
            Sig = sprintf('p_{%s - all}=%.3f',P_suffix,P);
        end
        t = text(-.025,-1,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
    end
    
    clear List_B
    
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
        
        %         for iSub = 1:size(B_rank,2)
        %             plot(...
        %                 mean(mean(B_rank(2:end,iSub,CdtToAvg,2,Eig_or_Mean),6),3),...
        %                 mean(Seed_Act_rank(iSub,2:end,Cdt_range_to_plot),4),...
        %                 'color',COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        %         end
        
        %         t=errorbar(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),mean(Seed_Act_rank(:,2:end)),....
        %             nansem(Seed_Act_rank(:,2:end))*0,nansem(Seed_Act_rank(:,2:end))*0,...
        %             nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1.5, 'markersize', 2,'color', 'k')
        
        plot(...
            mean(mean(mean(B_rank(2:end,:,CdtToAvg,2,Eig_or_Mean),6),3),2),...
            mean(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot),4)),'linewidth',1.5, ...
            'color', 'k')
        
        t=herrorbar(...
            mean(mean(mean(B_act(:,CdtToAvg,2,Eig_or_Mean),5),2)),...
            Y_act,...
            nansem(mean(mean(B_act(:,CdtToAvg,2,Eig_or_Mean),5),2)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        t=herrorbar(...
            mean(mean(mean(B_deact(:,CdtToAvg,2,Eig_or_Mean),5),2)),...
            Y_deact,...
            nansem(mean(mean(B_deact(:,CdtToAvg,2,Eig_or_Mean),5),2)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        t=herrorbar(...
            mean(mean(mean(B_all(:,CdtToAvg,2,Eig_or_Mean),5),2)),...
            Y_all,...
            nansem(mean(mean(B_all(:,CdtToAvg,2,Eig_or_Mean),5),2)),'.b');
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        List_B(:,iCdt) = mean(mean(B_all(:,CdtToAvg,2,Eig_or_Mean),5),2);
        
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
                Sig = 'p_{inter - act / deact}<0.001';
            else
                Sig = sprintf('p_{inter - act / deact}=%.3f',P);
            end
            t = text(.025,Y_act,sprintf(Sig));
            set(t,'fontsize',8);
            
            if P<0.05
                set(t,'fontsize',10,'fontweight','bold')
                set(t,'color','r');
            end
        end
        clear tmp
        
        % compute perm test Y_all
        tmp = mean(mean(B_all(:,CdtToAvg,2,Eig_or_Mean),5),2);
        
        for iPerm = 1:size(ToPermute,1)
            tmp3 = ToPermute(iPerm,:);
            tmp3 = repmat(tmp3',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp3);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = 'p_{all}<0.001';
        else
            Sig = sprintf('p_{all}=%.3f',P);
        end
        t = text(.025,Y_all,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp3
        
        if V1_2_A1==1
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -4 19.5])
        else
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-2:2:6 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-2';'0';'2';'4';'6'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -2 13.5])
        end
        
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
    
    for i=1:2
        switch i
            case 1
                subplot(3,3,6)
                P_suffix = 'Att_A > Att_V';
                tmp = List_B(:,1)-List_B(:,2);
            case 2
                subplot(3,3,8)
                P_suffix = 'Stim A > V';
                tmp = List_B(:,3)-List_B(:,4);
        end
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = sprintf('p_{%s - all}<0.001',P_suffix);
        else
            Sig = sprintf('p_{%s - all}=%.3f',P_suffix,P);
        end
        t = text(-.025,-1,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
    end
    
    clear List_B
    
    
    mtit(['Connectivity ' strrep(ROI_Suffix,'_2_',' to ') ' - ' Suffix ' - slope averaging'],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectivity_' ROI_Suffix '_CV_slope_averaging - ' Suffix ' .tif']), '-dtiff')
    
end

%% Plot slope for each of the 2 X 2 conditions and the slope resulting of their concatenation
% this is always plotted against the value for V-Fix averaged over
% attention (based on Cst S param at each vertex)

close all
for Eig_or_Mean = 1
    
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
        
        %         for iSub = 1:size(B_rank,2)
        %             plot(...
        %                 mean(B_rank(2:end,iSub,iCdt,2,Eig_or_Mean),6),...
        %                 mean(Seed_Act_rank(iSub,2:end,Cdt_range_to_plot),4),...
        %                 'color',COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        %         end
        
        %         t=errorbar(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),mean(Seed_Act_rank(:,2:end)),....
        %             nansem(Seed_Act_rank(:,2:end))*0,nansem(Seed_Act_rank(:,2:end))*0,...
        %             nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),nansem(B_rank(2:end,:,iCdt,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1, 'markersize', 2,'color', Colors(iCdt,:))
        
        plot(...
            mean(mean(B_rank(2:end,:,iCdt,2,Eig_or_Mean),6),2),...
            mean(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot),4)),...
            'linewidth',1.5, 'color', Colors(iCdt,:))
        
        t=herrorbar(...
            mean(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),...
            Y_act,...
            nansem(mean(B_act(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),...
            Y_deact,...
            nansem(mean(B_deact(:,iCdt,2,Eig_or_Mean),5)),'.b');
        set(t, 'linewidth', 1, 'markersize', 10, 'color', Colors(iCdt,:))
        t=herrorbar(...
            mean(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(B_all(:,iCdt,2,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        plot([-50 50],[0 0], '--k')
        plot([0 0],[-50 50], '--k')
        
        List_B(:,iCdt) = mean(B_all(:,iCdt,2,Eig_or_Mean),5);
        
        
        if V1_2_A1==1
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -4 19.5])
        else
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-2:2:6 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-2';'0';'2';'4';'6'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.05 .05 -2 13.5])
        end
        
        if iCdt==1 || iCdt==4
            if V1_2_A1
                ylabel('[V-Fix]')
            else
                ylabel('[A-Fix]')
            end
        end
        
        title(sprintf('Connectivity ; slope\n%s',Conditions_Names{iCdt}))
        
    end
    
    %     i. Att A > Att V
    %     ii.  Stim A > V
    %     iii. Stim A > V for AttA > Att V
    List_B(:,3) = [];
    
    for i=1:3
        switch i
            case 1
                subplot(3,3,2)
                P_suffix = 'Att_A > Att_V';
                tmp = sum(List_B(:,[1 2]),2)-sum(List_B(:,[3 4]),2);
            case 2
                subplot(3,3,4)
                P_suffix = 'Stim A > V';
                tmp = sum(List_B(:,[1 3]),2)-sum(List_B(:,[2 4]),2);
            case 3
                subplot(3,3,5)
                P_suffix = 'Stim A > V for Att_A > Att_V';
                tmp = diff(List_B(:,[1 3]),[],2)-diff(List_B(:,[2 4]),[],2);
        end
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = sprintf('p_{%s - all}<0.001',P_suffix);
        else
            Sig = sprintf('p_{%s - all}=%.3f',P_suffix,P);
        end
        t = text(-.025,-1,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
    end
    
    clear List_B
    
    
    for iCdt = 1:4
        
        switch iCdt
            case 1
                subplot(3,3,3)
                RankData = B_rank_att;
                ActData = B_act_att;
                DeactData = B_deact_att;
                AllData = B_all_att;
                CdtToPlot = 1;
            case 2
                subplot(3,3,6)
                RankData = B_rank_att;
                ActData = B_act_att;
                DeactData = B_deact_att;
                AllData = B_all_att;
                CdtToPlot = 2;
            case 3
                subplot(3,3,7)
                RankData = B_rank_sens;
                ActData = B_act_sens;
                DeactData = B_deact_sens;
                AllData = B_all_sens;
                CdtToPlot = 1;
            case 4
                subplot(3,3,8)
                RankData = B_rank_sens;
                ActData = B_act_sens;
                DeactData = B_deact_sens;
                AllData = B_all_sens;
                CdtToPlot = 2;
        end
        
        hold on
        grid on
        
        %         for iSub = 1:size(RankData,2)
        %             plot(...
        %                 mean(RankData(2:end,iSub,CdtToPlot,2,Eig_or_Mean),6),...
        %                 mean(Seed_Act_rank(iSub,2:end,Cdt_range_to_plot),4),...
        %                 'color',COLOR_Subject(iSub,:));%[.5 .5 .5]) % [.3 .3 1]
        %         end
        
        %         t=errorbar(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),mean(Seed_Act_rank(:,2:end)),....
        %             nansem(Seed_Act_rank(:,2:end))*0,nansem(Seed_Act_rank(:,2:end))*0,...
        %             nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),nansem(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),2),'o-b');
        %         set(t, 'linewidth', 1.5, 'markersize', 2,'color', 'k')
        
        plot(...
            mean(mean(RankData(2:end,:,CdtToPlot,2,Eig_or_Mean),6),2),...
            mean(mean(Seed_Act_rank(:,2:end,Cdt_range_to_plot),4)),'linewidth',1.5, ...
            'color', 'k')
        
        t=herrorbar(...
            mean(mean(ActData(:,CdtToPlot,2,Eig_or_Mean),5)),...
            Y_act,...
            nansem(mean(ActData(:,CdtToPlot,2,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        t=herrorbar(...
            mean(mean(DeactData(:,CdtToPlot,2,Eig_or_Mean),5)),...
            Y_deact,...
            nansem(mean(DeactData(:,CdtToPlot,2,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
        t=herrorbar(...
            mean(mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5)),'.b');
        
        switch iCdt
            case 1
                title(fprintf('Connectivity ; slope\nA attention'))
            case 2
                title(fprintf('Connectivity ; slope\nV attention'))
            case 3
                title(fprintf('Connectivity ; slope\nA stim'))
            case 4
                title(fprintf('Connectivity ; slope\nV stim'))
        end
        fprintf('\nmean +/- SEM: %f +/- %f\n\n', ...
            mean(mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5)),...
            nansem(mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5)));
        
        plot([-15 15],[0 0], '--k')
        plot([0 0],[-15 15], '--k')
        
        List_B(:,iCdt) = mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5);

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
                Sig = 'p_{inter - act / deact}<0.001';
            else
                Sig = sprintf('p_{inter - act / deact}=%.3f',P);
            end
            t = text(.025,Y_act,sprintf(Sig));
            set(t,'fontsize',8);
            
            if P<0.05
                set(t,'fontsize',10,'fontweight','bold')
                set(t,'color','r');
            end
        end
        clear tmp
        
        % compute perm test Y_all
        tmp = mean(AllData(:,CdtToPlot,2,Eig_or_Mean),5);
        
        for iPerm = 1:size(ToPermute,1)
            tmp3 = ToPermute(iPerm,:);
            tmp3 = repmat(tmp3',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp3);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = 'p_{all}<0.001';
        else
            Sig = sprintf('p_{all}=%.3f',P);
        end
        t = text(.025,Y_all,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp3
        
        if V1_2_A1==1
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-8:4:12 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-8';'-4';'0';'4';'8';'12'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.1 .05 -4 19.5])
        else
            set(gca,'ygrid','off','xtick',-3:.025:3,'xticklabel',-3:.025:3, ...
                'ytick', [-2:2:6 Y_deact Y_act Y_all],  ...
                'yticklabel',{'-2';'0';'2';'4';'6'; 'deactivated'; 'activated'; 'all vertices'})
            axis([-.135 .05 -2 13.5])
        end
        
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
    
    for i=1:2
        switch i
            case 1
                subplot(3,3,6)
                P_suffix = 'Att_A > Att_V';
                tmp = List_B(:,1)-List_B(:,2);
            case 2
                subplot(3,3,8)
                P_suffix = 'Stim A > V';
                tmp = List_B(:,3)-List_B(:,4);
        end
        for iPerm = 1:size(ToPermute,1)
            tmp2 = ToPermute(iPerm,:);
            tmp2 = repmat(tmp2',1,size(tmp,2));
            Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
        end
        P = sum( ...
            abs( Perms ) > ...
            repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
            / size(Perms,1);
        
        Sig = [];
        if P<0.001
            Sig = sprintf('p_{%s - all}<0.001',P_suffix);
        else
            Sig = sprintf('p_{%s - all}=%.3f',P_suffix,P);
        end
        t = text(-.025,-1,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
    end
    
    clear List_B
    
    mtit(['Connectivity ' strrep(ROI_Suffix,'_2_',' to ') ' - ' Suffix ' - concatenation'],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectivity_' ROI_Suffix '_CV_concatenation - ' Suffix ' .tif']), '-dtiff')
    
end
