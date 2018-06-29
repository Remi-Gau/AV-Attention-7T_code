function BOLD_vol_AllCdt_Plot
% Compare the different ROIs

clear; clc; close all;

Visibility='on';

NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

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

StartFolder = fullfile(pwd, '..', '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ResultsFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2/Figures/6_layers/All_cdtions';


ROI(1) = struct('name', 'A1_V_Act_surf', 'figname', 'A1: part activated by Visual stimuli');
ROI(end+1) = struct('name', 'PT_V_Act_surf', 'figname', 'PT: part activated by Visual stimuli');
ROI(end+1) = struct('name', 'V1_A_Act_surf', 'figname', 'V1: part activated by Audio stimuli');
ROI(end+1) = struct('name', 'V23_A_Act_surf', 'figname', 'V23: part activated by Audio stimuli');

ROI(end+1) = struct('name', 'A1_V_Deact_surf', 'figname', 'A1: part deactivated by Visual stimuli');
ROI(end+1) = struct('name', 'PT_V_Deact_surf', 'figname', 'PT: part deactivated by Visual stimuli');
ROI(end+1) = struct('name', 'V1_A_Deact_surf', 'figname', 'V1: part deactivated by Audio stimuli');
ROI(end+1) = struct('name', 'V23_A_Deact_surf', 'figname', 'V23: part deactivated by Audio stimuli');

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


for iROI=1:numel(ROI)
    
    LEGEND={sprintf('%s\nVolume analysis',ROI(iROI).figname)};
    Name =  ROI(iROI).name;
    
    ROIs_DATA.DATA_vol = squeeze(All_BOLD_Profile_vol(:,:,iROI,:));
    ROIs_DATA.DATA_vol_STD = squeeze(All_BOLD_Profile_vol_STD(:,:,iROI,:));
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

Transparent=1;

tmp = ROIs_DATA.DATA_vol-ROIs_DATA.DATA_vol_STD;
MIN = min(tmp(:));
tmp = ROIs_DATA.DATA_vol+ROIs_DATA.DATA_vol_STD;
MAX = max(tmp(:));


%%

figure('name', [Name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

iSubplot = 1;

for icdt=1:size(MEAN_vol,2)
    
    subplot(2,3,iSubplot)
    hold on
    grid on
    
    errorbar([1:size(MEAN_vol,1)]+.5, MEAN_vol(:,icdt), SEM_vol(:,icdt), ...
        'LineWidth', 2, 'Color', 'k')
    
    for iSubj=1:size(ROIs_DATA.DATA_vol,3)
        plot([1:size(MEAN_vol,1)]+.5, ROIs_DATA.DATA_vol(:,icdt,iSubj), '--','color', [.5 .5 .5])
    end
    
    plot([0 10],[0 0], '--k','LineWidth', 1)
    
    axis([0.9 7.1 MIN MAX])
    
    t=ylabel('Beta values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 1:7 ,'xticklabel', round(linspace(0,100,7)), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize-2)
    
    
    iSubplot=iSubplot+1;
end

subplot(2,3,1)
t=ylabel('A attention');
set(t,'fontsize',Fontsize+4);
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
set(t,'fontsize',Fontsize+4);

mtit(strrep(LEGEND{1},'_','-'), 'fontsize', 10, 'xoff',0,'yoff',.025)

print(gcf, fullfile(ResultsFolder,[Name, '_profiles_vol.tif']), '-dtiff')



end