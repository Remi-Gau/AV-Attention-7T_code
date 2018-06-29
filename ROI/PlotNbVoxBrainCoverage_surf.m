%% Compile values from different subjects
clear

StartFolder=fullfile(pwd, '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %%  '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %  '14';...
    '15';...
    '16'
    ];

AllSubjROI(1) = struct('name', [], 'size', [], 'Tab', []);
AllSubjROI(end+1) = struct('name', [], 'size', [], 'Tab', []);
AllSubjROI(end+1) = struct('name', [], 'size', [], 'Tab', []);
AllSubjROI(end+1) = struct('name', [], 'size', [], 'Tab', []);

for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Surfaces');
    
    load(fullfile(Results_Folder,['Subj_' SubjID '_brain_coverage.mat']),'Tab')
    load(fullfile(SubjectFolder,'ROI_MIPAV',['Subj_' SubjID '_ROI_coverage.mat']),'ROI')
    
    Coverage(:,:,SubjInd) = Tab(:,2:3);
    
    for iROI = 1:numel(AllSubjROI)
        AllSubjROI(iROI).name = ROI(iROI).name;
        AllSubjROI(iROI).Tab(:,:,SubjInd) = ROI(iROI).Tab(:,2:3);
        AllSubjROI(iROI).size(SubjInd,1:2) = cellfun('length',(ROI(iROI).VertOfInt));
    end
    
end


%% Figure Nb Voxels

close all

COLOR_Subject= [
    127,127,127;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

FigDim = [100 100 1500 1000];
Visible = 'on';

Scatter = linspace(1.1,1.5,size(SubjectList,1));


%%
figure('position', FigDim, 'name', 'ROI size', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

iSubPlot = 1;

for iROI = 1:numel(AllSubjROI)
    
    subplot(numel(AllSubjROI),2,iSubPlot)
    hold on
    errorbar(1,mean(AllSubjROI(iROI).size(:,1)), nansem(AllSubjROI(iROI).size(:,1)), '.k')
    for iSubj=1:size(SubjectList,1)
        plot(Scatter(iSubj),AllSubjROI(iROI).size(iSubj,1), 'marker', '.', 'markersize', 30, ....
            'color', COLOR_Subject(iSubj,:))
    end
    axis([0.9 1.6 0 max(AllSubjROI(iROI).size(:))])
    set(gca,'xtick',Scatter, 'xticklabel', SubjectList, 'ygrid', 'on')
    ylabel(strrep(AllSubjROI(iROI).name, '_', ' '))
    
    iSubPlot = iSubPlot + 1;
    
    subplot(numel(AllSubjROI),2,iSubPlot)
    hold on
    errorbar(1,mean(AllSubjROI(iROI).size(:,2)), nansem(AllSubjROI(iROI).size(:,2)), '.k')
    for iSubj=1:size(SubjectList,1)
        plot(Scatter(iSubj),AllSubjROI(iROI).size(iSubj,2), 'marker', '.', 'markersize', 30, ....
            'color', COLOR_Subject(iSubj,:))
    end
    axis([0.9 1.6 0 max(AllSubjROI(iROI).size(:))])
    set(gca,'xtick',Scatter, 'xticklabel', SubjectList, 'ygrid', 'on')
    
    iSubPlot = iSubPlot + 1;
    
end

subplot(numel(AllSubjROI),2,1)
title('LEFT')
subplot(numel(AllSubjROI),2,2)
title('RIGHT')

mtit('ROI size', 'xoff',0,'yoff',.025)

print(gcf, fullfile(StartFolder, 'Figures', 'Profiles', ...
    'NbVoxelLeftRightROI_surf.tif'), '-dtiff')


%%
FigName={'Percent ROI not covered', 'Percent ROI small values', 'Percent ROI with values'};

for iTab=1:3
    
    figure('position', FigDim, 'name', FigName{iTab}, 'Color', [1 1 1], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    iSubPlot = 1;
    
    for iROI = 1:numel(AllSubjROI)
        
        subplot(numel(AllSubjROI),2,iSubPlot)
        hold on
        errorbar(1,mean(AllSubjROI(iROI).Tab(iTab,1,:)), nansem(AllSubjROI(iROI).Tab(iTab,1,:)), '.k')
        for iSubj=1:size(SubjectList,1)
            plot(Scatter(iSubj),AllSubjROI(iROI).Tab(iTab,1,iSubj), 'marker', '.', 'markersize', 30, ....
                'color', COLOR_Subject(iSubj,:))
        end
        axis([0.9 1.6 0 1])
        set(gca,'xtick',Scatter, 'xticklabel', SubjectList, 'ygrid', 'on')
        ylabel(strrep(AllSubjROI(iROI).name, '_', ' '))
        
        iSubPlot = iSubPlot + 1;
        
        subplot(numel(AllSubjROI),2,iSubPlot)
        hold on
        errorbar(1,mean(AllSubjROI(iROI).Tab(iTab,2,:)), nansem(AllSubjROI(iROI).Tab(iTab,2,:)), '.k')
        for iSubj=1:size(SubjectList,1)
            plot(Scatter(iSubj),AllSubjROI(iROI).Tab(iTab,2,iSubj), 'marker', '.', 'markersize', 30, ....
                'color', COLOR_Subject(iSubj,:))
        end
        axis([0.9 1.6 0 1])
        set(gca,'xtick',Scatter, 'xticklabel', SubjectList, 'ygrid', 'on')
        
        iSubPlot = iSubPlot + 1;
        
    end
    
    subplot(numel(AllSubjROI),2,1)
    title('LEFT')
    subplot(numel(AllSubjROI),2,2)
    title('RIGHT')
    
    mtit(FigName{iTab}, 'xoff',0,'yoff',.025)
    
        print(gcf, fullfile(StartFolder, 'Figures', 'Profiles', ...
    [strrep(FigName{iTab},' ','_') '_surf.tif']), '-dtiff')
    
end