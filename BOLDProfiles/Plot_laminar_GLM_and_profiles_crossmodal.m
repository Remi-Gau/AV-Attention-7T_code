clc; clear; close all;

PlotSubjects = 0;
Median = 1;
Switch = 1;
NbLayers = 6;

StartDirectory=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartDirectory, 'SubFun')))
Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

SourceFolder = fullfile(StartDirectory, 'Figures', 'ProfilesSurface', strcat(num2str(NbLayers), '_layers'));

FigureFolder = fullfile(StartDirectory, 'Figures', strcat(num2str(NbLayers+2), '_layers'));
mkdir(FigureFolder)

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '15';...
    '16'
    ];

ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3';...
    'V1_act';...
    'V1_deact';...
    'V23_act';...
    'V23_deact';...
    };

FigDim = [100 100 900 500];
Visible = 'on';
Transparent = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);


for iSubj=1:size(SubjectList,1)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


if Median
    MedianSufix = '_Median'; %#ok<*UNRCH>
else
    MedianSufix = '';
end

if PlotSubjects
    SubjSuffix = '_Subj'; %#ok<*UNRCH>
else
    SubjSuffix = '';
end



%% Get data for BOLD
load(fullfile(SourceFolder, strcat('Data_Surf', MedianSufix ,'_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% 1 A against baseline under A att
% 2 V against baseline under A att
% 3 AV against baseline under A att
% 4 A against baseline under V att
% 5 V against baseline under V att
% 6 AV against baseline under V att
Target=1;
for iCond = [1:3 5:7]
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).MainEffects.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% 7 AV-A under A att
% 8 AV-A under V att
% 9 AV-A under A-V att
% 10 AV-V under A att
% 11 AV-V under V att
% 12 AV-V under A-V att
Target=7;
for iCond = 7:12
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).BiVSUniSep.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end


Include =[];
for iROI=1:numel(ROIs)
    if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
        temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
        Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
    end
end
clear AllSubjects_Data temp iROI


%% Plot AV-A and AV-V

close all
clear DATA


DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 1;
DATA.PlotInset = 0;


for iROI = 5:8
    
    figure('position', FigDim, 'name', 'AV-A & AV-V', ...
        'Color', [1 1 1], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    DATA.MVPA = 0;
    DATA.YLabelInset = 1;
    DATA.OneSideTTest = {'both' 'both' 'both'};
    
    SubPlotCurves = [1:3 7:9];
    SubPlotLamGLM = [4:6 10:12];
    COI = 7:12;
    
    for Cdt2Plot = 1:6
        
        if iROI==6 || iROI==8
            AxisLim = 3;
        elseif iROI<3 && Cdt2Plot<4 || iROI>2 && Cdt2Plot>3
            AxisLim = 1;
        else
            AxisLim = 2;
        end
        
        if AxisLim==1
            DATA.InsetLim = [1 .4;-1 -.4];
            if DATA.WithSubj
                DATA.MIN = -1.5;
                DATA.MAX = 1.5;
            else
                DATA.MIN = -.5;
                DATA.MAX = .5;
            end
        elseif AxisLim==2
            DATA.InsetLim = [5 1;-2 -0.4];
            if DATA.WithSubj
                DATA.MIN = -1.5;
                DATA.MAX = 4;
            else
                DATA.MIN = -1;
                DATA.MAX = 5;
            end
        elseif AxisLim==3
            DATA.InsetLim = [2.4 .5;-2.4 -.5];
            if DATA.WithSubj
                DATA.MIN = -1.5;
                DATA.MAX = 4;
            else
                DATA.MIN = -1.2;
                DATA.MAX = 1.2;
            end
        end
        
        iCond = COI(Cdt2Plot);
        
        subplot(4,3,SubPlotCurves(Cdt2Plot))
        PlotRectangle(NbLayers,Fontsize,Switch)
        subplot(4,3,SubPlotCurves(Cdt2Plot))
        
        DATA.Name = '';
        DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
        DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
        DATA.Color =  [0 0 0];
        DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
        
        PlotProfileAndBetas(DATA)
        
        ax = subplot(4,3,SubPlotLamGLM(Cdt2Plot));
        axis('off')
        DATA.ax = ax.Position;
        DATA.ToPermute = ToPermute;
        PlotInsetFinal(DATA)
        
    end
    
    subplot(4,3,1)
    title('Attention A')
    ylabel('AV-A')
    subplot(4,3,2)
    title('Attention V')
    subplot(4,3,3)
    title('Attention A - Attention V')
    subplot(4,3,7)
    ylabel('AV-V')
    
    mtit(ROIs{iROI},'xoff', 0, 'yoff', +0.04, 'fontsize', 14)
    
    print(gcf, fullfile(FigureFolder,['Crossmodal_effects_' ROIs{iROI} '_' MedianSufix SubjSuffix '.tif']), '-dtiff')
    
end

return


%% Plot AV-A and AV-V all ROIs in one
close all
clear DATA

FigDim = [100 100 1800 1000];

DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 1;
DATA.PlotInset = 0;

figure('position', FigDim, 'name', 'AV-A & AV-V', ...
    'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.YLabelInset = 1;
DATA.OneSideTTest = {'both' 'both' 'both'};

Subplot_first = 1:6:24;

for iROI = 1:4
    
    if iROI<3
        COI = 7:9;
    else
        COI = 10:12;
    end
    
    Subplot = Subplot_first(iROI);
    
    for iCond = COI
        
        DATA.InsetLim = [1.8 .4;-1.4 -.4];
        if DATA.WithSubj
            DATA.MIN = -1.5;
            DATA.MAX = 2;
        else
            DATA.MIN = -.5;
            DATA.MAX = .7;
        end
        
        DATA.Name = '';
        DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
        DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
        DATA.Color =  [0 0 0];
        DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
        
        
        Subplot = Subplot + 3;
        ax = subplot(8,3,Subplot);
        axis('off')
        DATA.ax = ax.Position;
        DATA.ax(2) = DATA.ax(2)-0.01;
        DATA.ax(4) = DATA.ax(4)-0.02;
        DATA.ToPermute = ToPermute;
        PlotInsetFinal(DATA)

        Subplot = Subplot - 3;
        subplot(8,3,Subplot)
        PlotRectangle(NbLayers,Fontsize,Switch)
        subplot(8,3,Subplot)

        PlotProfileAndBetas(DATA)
        set(gca,'fontsize', 11)

        Subplot = Subplot + 1;
        
    end
    
    
%     subplot(8,3,1)
%     t=ylabel(sprintf('A1\nParam. est. [a u]'))
%     set(t, 'fontsize', 10) 
%     subplot(8,3,7)
%     t=ylabel(sprintf('PT\nParam. est. [a u]'))
%     set(t, 'fontsize', 10) 
%     subplot(8,3,13)
%     t=ylabel(sprintf('V1\nParam. est. [a u]'))
%     set(t, 'fontsize', 10) 
%     subplot(8,3,19)
%     t=ylabel(sprintf('V2/3\nParam. est. [a u]'))
%     set(t, 'fontsize', 10) 
    
    
%     subplot(8,3,4)
%     t=ylabel('Param. est. [a u]')
%     set(t, 'fontsize', 8)
%     subplot(8,3,10)
%     t=ylabel('Param. est. [a u]')
%     set(t, 'fontsize', 8)
%     subplot(8,3,16)
%     t=ylabel('Param. est. [a u]')
%     set(t, 'fontsize', 8)
%     subplot(8,3,22)
%     t=ylabel('Param. est. [a u]')
%     set(t, 'fontsize', 8)
    
    
%     subplot(8,3,1)
%     t=title('(AV-A)_{A att}')
%     set(t, 'fontsize', 10)
%     t=subplot(8,3,2)
%     title('(AV-A)_{V att}')
%     set(t, 'fontsize', 10)    
%     t=subplot(8,3,3)
%     title('(AV-A)_{A att - V att}')
%     set(t, 'fontsize', 10)    
%     
%     subplot(8,3,13)
%     t=title('(AV-V)_{A att}')
%     set(t, 'fontsize', 10)
%     subplot(8,3,14)
%     t=title('(AV-V)_{V att}')
%     set(t, 'fontsize', 10)
%     subplot(8,3,15)
%     t=title('(AV-V)_{A att - V att}')
%     set(t, 'fontsize', 10)
    
    %     mtit(ROIs{iROI},'xoff', 0, 'yoff', +0.04, 'fontsize', 14)
    
        print(gcf, fullfile(FigureFolder,['Crossmodal_effects_All_ROIs_' MedianSufix SubjSuffix '.tif']), '-dtiff')
    
end





%% All conditions
close all
clear DATA

Transparent = 1;
Fontsize = 13;

figure('position', FigDim, 'name', 'AllCdtions', ...
    'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);


DATA.YLabelInset = 1;
DATA.OneSideTTest = {'both' 'both' 'both'};

Subplot = 1;

for iROI = 1:4
    
    for iCond = 1:3
        
        if iROI<3 && iCond==2
            AxisLim = 1;
        elseif iROI>2 && iCond==1
            AxisLim = 2;
        elseif iROI<3
            AxisLim = 3;
        else
            AxisLim = 4;
        end
        
        if AxisLim==1
            MIN = -.4;
            MAX = .2;
            if PlotSubjects
                MIN = -.8;
                MAX = .8;
            end
        elseif AxisLim==2
            MIN = -1.6;
            MAX = .05;
            if PlotSubjects
                MIN = -2;
                MAX = .3;
            end
        elseif AxisLim==3
            MIN = -.1;
            MAX = 3.2;
            if PlotSubjects
                MIN = -.3;
                MAX = 3;
            end
        elseif AxisLim==4
            MIN = -.1;
            MAX = 2.5;
            if PlotSubjects
                MIN = -.3;
                MAX = 3;
            end
        end
        
        
        subplot(4,3,Subplot)
        if iROI==4
            PlotRectangle(NbLayers,Fontsize,Switch)
            subplot(4,3,Subplot)
        end
        
        % box on
        hold on
        grid on
        
        plot([0 10], [0 0], ':k', 'LineWidth', 1.5)
        
        Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
        if PlotSubjects
            for SubjInd = 1:size(Data,1)
                plot(1:NbLayers, Data(SubjInd,:), '-', ...
                    'LineWidth', .5, 'Color', [27,158,119]/256+.3);
            end
        end
        
        Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond+3,iROI,:)));
        if PlotSubjects
            for SubjInd = 1:size(Data,1)
                plot(1:NbLayers, Data(SubjInd,:), '-', ...
                    'LineWidth', .5, 'Color', [117,112,179]/256+.3);
            end
        end
        
        Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
        shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
            {'Marker', '.', 'MarkerSize', 12, 'LineWidth', 2, 'Color', 'k'}, Transparent)
        
        Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond+3,iROI,:)));
        shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
            {'Marker', '.', 'MarkerSize', 12, 'linestyle', ':', 'LineWidth', 2, 'Color', 'k'}, Transparent)
        
        axis([0.5 NbLayers+.5 MIN MAX])
        
        if iCond==1
            %             ylabel(sprintf('%s\nParam. est. [a u]', ROIs{iROI}))
        end
        
        Subplot = Subplot +1;
        
        set(gca, 'xtick', [])
        
    end
    
    
    
end

% subplot(4,3,1)
% title('A')
% subplot(4,3,2)
% title('V')
% subplot(4,3,3)
% title('AV')

print(gcf, fullfile(FigureFolder,['All_effects_' MedianSufix SubjSuffix '.tif']), '-dtiff')


