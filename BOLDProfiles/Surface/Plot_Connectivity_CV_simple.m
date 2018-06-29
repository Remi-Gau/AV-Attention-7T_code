%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

COLOR_Subject = ColorSubject();

V1_2_A1 = 1;
if V1_2_A1==1
    ROI_Suffix = 'V1_2_A1';
else
    ROI_Suffix = 'A1_2_V1';
end
Y_all = 1;

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

load(fullfile(Results_Folder,strcat('Connectivity_all_', ROI_Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')))

B_all_V1toA1 = B_all;

[B_all_V1toA1(:,1,1,2) 1./B_all(:,1,1,2)]

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
        
        t=herrorbar(...
            mean(mean(B_all(:,iCdt,1,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(B_all(:,iCdt,1,Eig_or_Mean),5)),'.k');
        set(t, 'linewidth', 1, 'markersize', 10)
                plot([-6 6],[0 0], '--k')
                plot([0 0],[-2 2], '--k')
        
        List_B(:,iCdt) = mean(B_all(:,iCdt,1,Eig_or_Mean),5);
        
        
        set(gca,'ygrid','off','xtick',-6:.5:6,'xticklabel',-6:.5:6, ...
            'ytick', -6:1:6,  ...
            'yticklabel',-6:1:6)
        axis([-1.25 1.25 -2 2])
        
        if iCdt==1 || iCdt==4
            if V1_2_A1
                ylabel('[V-Fix]')
            else
                ylabel('[A-Fix]')
            end
        end
        
        % compute perm test Y_all
        tmp = mean(B_all(:,iCdt,1,Eig_or_Mean),5);
        
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
        t = text(.75,1.5,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp3
        
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
        t = text(-.75,-1,sprintf(Sig));
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
                AllData = B_all_att;
                CdtToPlot = 1;
            case 2
                subplot(3,3,6)
                AllData = B_all_att;
                CdtToPlot = 2;
            case 3
                subplot(3,3,7)
                AllData = B_all_sens;
                CdtToPlot = 1;
            case 4
                subplot(3,3,8)
                AllData = B_all_sens;
                CdtToPlot = 2;
        end
        
        hold on
        grid on
        
        t=herrorbar(...
            mean(mean(AllData(:,CdtToPlot,1,Eig_or_Mean),5)),...
            Y_all,...
            nansem(mean(AllData(:,CdtToPlot,1,Eig_or_Mean),5)),'.b');
                plot([-6 6],[0 0], '--k')
                plot([0 0],[-2 2], '--k')
        
        List_B(:,iCdt) = mean(AllData(:,CdtToPlot,1,Eig_or_Mean),5);
        
        % compute perm test Y_all
        tmp = mean(AllData(:,CdtToPlot,1,Eig_or_Mean),5);
        
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
        t = text(.75,1.5,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
        clear tmp tmp3
        
        
        set(gca,'ygrid','off','xtick',-6:.5:6,'xticklabel',-6:.5:6, ...
            'ytick', -6:1:6,  ...
            'yticklabel',-6:1:6)
        axis([-1.25 1.25 -2 2])
        
        
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
        t = text(-.75,-1,sprintf(Sig));
        set(t,'fontsize',8);
        
        if P<0.05
            set(t,'fontsize',10,'fontweight','bold')
            set(t,'color','r');
        end
    end
    
    clear List_B
    
    mtit(['Connectivity ' strrep(ROI_Suffix,'_2_',' to ') ' - ' Suffix ' - concatenation'],'xoff',0,'yoff',.04,'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Connectivity_' ROI_Suffix '_concatenation - ' Suffix ' .tif']), '-dtiff')
    
end
