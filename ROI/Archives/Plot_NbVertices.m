clear; close all; clc;

NbLayers = 7;
NbLayers = NbLayers+1;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

StartFolder=pwd;
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ROIs_Ori = {...
%     'A1';...
%     'TE1.0';...
%     'TE1.1';...
%     'TE1.2';...
%     'STG_Post';...
    'V1v';...
    'V1d';...
    'V2v';...
    'V2d';...
    'V3v';...
    'V3d';...
    'V4';...
    'V5';...
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
    Results_Folder = fullfile(SubjectFolder, 'Transfer', 'ProfilesSurface', TargetLayerFile);
    
    File2Load = fullfile(Results_Folder, 'ROIs.mat');
    
    if exist(File2Load, 'file')
        
        load(File2Load, 'ROIs')
        
        for iROI=1:size(ROIs_Ori,1)

            id = ismember(ROIs(:,1),ROIs_Ori(iROI));
            
            NbVertices(SubjInd,iROI) = numel(ROIs{id,3})+numel(ROIs{id,4}); %#ok<*SAGROW>
            
        end
        
    else
        
        NbVertices(SubjInd,:) = nan(numel(ROIs_Ori),1);
        
    end
end

%%
figure('name', 'Number of vertices per ROI', 'position', [100 100 1500 1000])

hold on

for iROI = 1:numel(ROIs_Ori)
    
    for iSubj=1:size(NbVertices,1)
        plot(iROI+Jitter(iSubj), NbVertices(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-.1, nanmean(NbVertices(:,iROI)), nansem(NbVertices(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

tmp = char(ROIs_Ori);
% tmp = strrep(tmp, '_', ' ')

% axis([0.5 size(tmp,1)+.5 -100 23000])

set(gca,'tickdir', 'out', ...
    'xtick', 1:(size(tmp,1)), ...
    'xticklabel', tmp,...
    'fontsize', Fontsize)

t = xlabel('ROIs');
set(t, 'fontsize', Fontsize)
t = ylabel('Nb vertices');
set(t, 'fontsize', Fontsize)

print(gcf, fullfile(FigureFolder,'NumberVerticesV.tif'), '-dtiff')
