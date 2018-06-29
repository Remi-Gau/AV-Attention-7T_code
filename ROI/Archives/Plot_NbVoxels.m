clear; close all; clc;

NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

StartFolder=pwd;
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ROIs = {...
%     'A1';...
%     'TE1.0';...
%     'TE1.1';...
%     'TE1.2';...
    'STG_Post';...
%     'V1v';...
%     'V1d';...
%     'V2v';...
%     'V2d';...
%     'V3v';...
%     'V3d';...
%     'V4';...
%     'V5';...
    };

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


FigureFolder = fullfile(StartFolder, 'Figures');

NbSubject = size(SubjectList,1);

Fontsize = 10;

Jitter = linspace(0.015,0.55,NbSubject);
Jitter = Jitter(randperm(NbSubject));


%% Gets data for each subject
for SubjInd = 1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    AnalysisFolder = fullfile(SubjectFolder, 'Transfer', 'Profiles', TargetLayerFile);
    
    for iROI=1:size(ROIs,1)
        
        File2Load = fullfile(AnalysisFolder, strcat('Data_Block_', ...
            ROIs{iROI,1}, '_', TargetLayerFile, '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'Data_ROI')
            NbVoxels(SubjInd,iROI) = Data_ROI.ROI_Coverage(2); %#ok<*SAGROW>
            Coverage(SubjInd,iROI) = Data_ROI.ROI_Coverage(1)/Data_ROI.ROI_Coverage(2);
        else
            NbVoxels(SubjInd,iROI) = NaN;
            Coverage(SubjInd,iROI) = NaN;
        end
        
    end
    
end

%%
figure('name', 'Number of vertices per ROI', 'position', [100 100 1500 1000])

hold on

for iROI = 1:numel(ROIs)
    
    for iSubj=1:size(NbVoxels,1)
        plot(iROI+Jitter(iSubj), NbVoxels(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-.1, nanmean(NbVoxels(:,iROI)), nansem(NbVoxels(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

tmp = char(ROIs);
tmp = strrep(tmp, '_', ' ')

% axis([0.5 size(tmp,1)+.5 -100 23000])

set(gca,'tickdir', 'out', ...
    'xtick', 1:(size(tmp,1)), ...
    'xticklabel', tmp,...
    'fontsize', Fontsize)

t = xlabel('ROIs');
set(t, 'fontsize', Fontsize)
t = ylabel('Nb voxels');
set(t, 'fontsize', Fontsize)

print(gcf, fullfile(FigureFolder,'NumberVoxelsSTG.tif'), '-dtiff')


%%
figure('name', 'Coverage', 'position', [100 100 1500 1000])

hold on

for iROI = 1:numel(ROIs)
    
    for iSubj=1:size(Coverage,1)
        plot(iROI+Jitter(iSubj), Coverage(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-.1, nanmean(Coverage(:,iROI)), nansem(Coverage(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

tmp = char(ROIs);
tmp = strrep(tmp, '_', ' ')

% axis([0.5 size(tmp,1)+.5 -100 23000])

set(gca,'tickdir', 'out', ...
    'xtick', 1:(size(tmp,1)), ...
    'xticklabel', tmp,...
    'fontsize', Fontsize)

t = xlabel('ROIs');
set(t, 'fontsize', Fontsize)
t = ylabel('Coverage');
set(t, 'fontsize', Fontsize)

print(gcf, fullfile(FigureFolder,'CoverageSTG.tif'), '-dtiff')

cd(StartFolder)





