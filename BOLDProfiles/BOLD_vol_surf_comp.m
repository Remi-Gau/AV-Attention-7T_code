function BOLD_vol_surf_comp
% Compare the different ROIs

clear; clc; close all;

Visibility='on';

NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '14';...
    '15';...
    '16'
    ];

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ResultsFolder = fullfile(StartFolder,'Figures','Profiles');


%%
ROI(1) = struct('name', 'A1_surf');
ROI(end+1) = struct('name', 'PT_surf_thres');

ROI(end+1) = struct('name', 'V1_surf_thres');
ROI(end+1) = struct('name', 'V2-3_surf_thres');

All_BOLD_Profile_vol = nan(NbLayers,6,numel(ROI),size(SubjectList,1));
All_BOLD_Profile_vol_STD = nan(NbLayers,6,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Volumes', TargetLayerFile);
    
    for ROI_Ind = 1:numel(ROI)
        File2Load = fullfile(SubjectFolder, ...
            strcat('Data_Block_', ROI(ROI_Ind).name, '_', TargetLayerFile, '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'Data_ROI')
            All_BOLD_Profile_vol(:,:,ROI_Ind,SubjInd) = Data_ROI.MEAN(2:NbLayers+1,:);
            All_BOLD_Profile_vol_STD(:,:,ROI_Ind,SubjInd) = squeeze(nanmean(Data_ROI.LayerSEM(2:NbLayers+1,:,:),2));
        else
            error('File %s is missing.', File2Load)
            All_BOLD_Profile_vol(:,:,ROI_Ind,SubjInd) = nan(NbLayers,6);  %#ok<*UNRCH>
            All_BOLD_Profile_vol_STD = nan(NbLayers,6);
        end
        
    end
end


%%
clear ROI
ROI(1) = struct('name', 'A1');
ROI(end+1) = struct('name', 'PT');

ROI(end+1) = struct('name', 'V1');
ROI(end+1) = struct('name', 'V2-3');

All_BOLD_Profile_surf = nan(NbLayers+2,6,numel(ROI),size(SubjectList,1));
All_BOLD_Profile_surf_STD = nan(NbLayers+2,6,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Surfaces');
    for ROI_Ind = 1:numel(ROI)
        File2Load = fullfile(SubjectFolder,...
            strcat('Data_Surf_Block_', ROI(ROI_Ind).name, '_', num2str(NbLayers+2), '_layers.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'Data_ROI')
            All_BOLD_Profile_surf(:,:,ROI_Ind,SubjInd) =  Data_ROI.MEAN;
            All_BOLD_Profile_surf_STD(:,:,ROI_Ind,SubjInd) = squeeze(nanmean(Data_ROI.LayerSEM,2));
        else
            error('File %s is missing.', File2Load)
            All_BOLD_Profile_surf(:,:,ROI_Ind,SubjInd) = nan(NbLayers,6);
            All_BOLD_Profile_surf_STD(:,:,ROI_Ind,SubjInd) = nan(NbLayers,6);
        end
        
    end
end


for iROI=1:numel(ROI)
    
    LEGEND={ROI(iROI).name}';
    Name =  ROI(iROI).name;
    
    ROIs_DATA.DATA_vol = squeeze(All_BOLD_Profile_vol(:,:,iROI,:));
    ROIs_DATA.DATA_vol_STD = squeeze(All_BOLD_Profile_vol_STD(:,:,iROI,:));
    ROIs_DATA.DATA_surf = squeeze(All_BOLD_Profile_surf(:,:,iROI,:));
    ROIs_DATA.DATA_surf_STD = squeeze(All_BOLD_Profile_surf_STD(:,:,iROI,:));
    ROIs_DATA.name = ROI(iROI).name;
    
    PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, Visibility, ResultsFolder)
    
    close all
end


end



function PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, Visibility, ResultsFolder)

Fontsize = 8;

% Color for Subjects
COLOR_Subject= [
    0,0,0;
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

mn = size(ROIs_DATA.DATA_vol,3);

MEAN_vol=nanmean(ROIs_DATA.DATA_vol,3);
STD_vol=nanstd(ROIs_DATA.DATA_vol,3);
SEM_vol=nanstd(ROIs_DATA.DATA_vol,3);

MEAN_surf=nanmean(ROIs_DATA.DATA_surf,3);
STD_surf=nanstd(ROIs_DATA.DATA_surf,3);
SEM_surf=nanstd(ROIs_DATA.DATA_surf,3);

Transparent=1;

tmp = ROIs_DATA.DATA_vol-ROIs_DATA.DATA_vol_STD;
tmp2 = ROIs_DATA.DATA_surf-ROIs_DATA.DATA_surf_STD;
MIN = min([tmp(:);tmp2(:)]);
tmp = ROIs_DATA.DATA_vol+ROIs_DATA.DATA_vol_STD;
tmp2 = ROIs_DATA.DATA_surf+ROIs_DATA.DATA_surf_STD;
MAX = max([tmp(:);tmp2(:)]);




%%

figure('name', [Name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

iSubplot = 1;

for icdt=1:size(MEAN_vol,2)
    
    subplot(2,3,iSubplot)
    hold on
    grid on
    
    shadedErrorBar(1:size(MEAN_vol,1), MEAN_vol(:,icdt), SEM_vol(:,icdt), ...
        {'LineWidth', 1, 'Color', 'k'}, Transparent)
    shadedErrorBar(0:size(MEAN_surf,1)-1, flipud(MEAN_surf(:,icdt)), flipud(SEM_surf(:,icdt)), ...
        {'LineStyle', '--', 'LineWidth', 1, 'Color', 'k'}, Transparent)
    
    plot([0 size(ROIs_DATA.DATA_surf,1)],[0 0], ':k','LineWidth', 1)
    
    axis([-0.1 size(ROIs_DATA.DATA_surf,1)-.9 MIN MAX])
        
    t=ylabel('Beta values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 0:size(ROIs_DATA.DATA_surf,1)-1 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_surf,1)))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize-2)
    
    iSubplot=iSubplot+1;
end

mtit(LEGEND{1})

print(gcf, fullfile(ResultsFolder, Name, [LEGEND{1}, '_comp_vol_surf.tif']), '-dtiff')


%%

for SubjInd=1:mn
    
    figure('name', [Name '-SUBJECTS-' SubjectList(SubjInd)], 'Position', [100, 100, 1500, 1000], ...
        'Color', [1 1 1], 'Visible', Visibility);
    
    iSubplot = 1;
    
%     tmp = ROIs_DATA.DATA_vol(:,:,SubjInd)-ROIs_DATA.DATA_vol_STD(:,:,SubjInd);
%     tmp2 = ROIs_DATA.DATA_surf(:,:,SubjInd)-ROIs_DATA.DATA_surf_STD(:,:,SubjInd);
%     MIN = min([tmp(:);tmp2(:)]);
%     tmp = ROIs_DATA.DATA_vol(:,:,SubjInd)+ROIs_DATA.DATA_vol_STD(:,:,SubjInd);
%     tmp2 = ROIs_DATA.DATA_surf(:,:,SubjInd)+ROIs_DATA.DATA_surf_STD(:,:,SubjInd);
%     MAX = max([tmp(:);tmp2(:)]);
    
    for icdt=1:size(MEAN_vol,2)
        
        subplot(2,3,iSubplot)
        hold on
        grid on
        
        shadedErrorBar(1:size(MEAN_vol,1), ROIs_DATA.DATA_vol(:,icdt,SubjInd), ...
            ROIs_DATA.DATA_vol_STD(:,icdt,SubjInd), ...
            {'LineWidth', 2 'Color', COLOR_Subject(SubjInd,:)}, Transparent)
        shadedErrorBar(0:size(MEAN_surf,1)-1, flipud(ROIs_DATA.DATA_surf(:,icdt,SubjInd)), ...
            flipud(ROIs_DATA.DATA_surf_STD(:,icdt,SubjInd)), ...
            {'LineStyle', '--', 'LineWidth', 1, 'Color', COLOR_Subject(SubjInd,:)}, Transparent)
        
        plot([0.5 size(ROIs_DATA.DATA_vol,1)+.5],[0 0], ':k','LineWidth',1)
        
        t=ylabel('Beta values');
        set(t,'fontsize',Fontsize); clear t
        
        t=xlabel('% Cortical depth');
        set(t,'fontsize',Fontsize); clear t
        
        set(gca,'tickdir', 'out', 'xtick', 0:size(ROIs_DATA.DATA_surf,1)-1 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_surf,1)))), ...
            'ticklength', [0.01 0.01], 'fontsize', Fontsize-2)
        
        axis([-0.1 size(ROIs_DATA.DATA_surf,1)-.9 MIN MAX])
        
        
        iSubplot=iSubplot+1;
    end
    
    subplot(2,3,1)
    t=ylabel('A attention');
    set(t,'fontsize',Fontsize+2);
    t=title('A stimulation');
    set(t,'fontsize',Fontsize);
    
    subplot(2,3,2)
    t=title('V stimulation');
    set(t,'fontsize',Fontsize);
    
    subplot(2,3,3)
    t=title('AV stimulation');
    set(t,'fontsize',Fontsize);
    
    subplot(2,3,4)
    t=ylabel('V attention');
    set(t,'fontsize',Fontsize+2);
    
    mtit(['Subject ' SubjectList(SubjInd,:) ' - ' LEGEND{1}], 'xoff', 0,'yoff',.025);
    
    
    print(gcf, fullfile(ResultsFolder, Name, ['Subject ' SubjectList(SubjInd,:) ' - ' LEGEND{1}, '_comp_vol_surf.tif']), '-dtiff')
end


end