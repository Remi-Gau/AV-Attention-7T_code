%%
clear; clc;

StartFolder = pwd;

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

ROIs = { ...
    'A_lh'; ...
    'A_rh'; ...
    'V_lh'; ...
    'V_rh'; ...
    };

ListSubROI = {...
    'A1_L'
    'TE_1.0_L_MNI'
    'TE_1.1_L_MNI'
    'TE_1.2_L_MNI'
    'PT_L_MNI'    % 5
    'V1_L_MNI'
    'V2_L_MNI'
    'V3d_L_MNI'
    'V3v_L_MNI'
    'V4d_L_MNI'
    'V4v_L_MNI'
    'V5_L_MNI'  % 12
    'A1_R'      % 13
    'TE_1.0_R_MNI'
    'TE_1.1_R_MNI'
    'TE_1.2_R_MNI'
    'PT_R_MNI'  % 17
    'V1_R_MNI'
    'V2_R_MNI'
    'V3d_R_MNI'
    'V3v_R_MNI'
    'V4d_R_MNI'
    'V4v_R_MNI'
    'V5_R_MNI'
    };

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

FigureFolder = fullfile(StartFolder, 'ROI_Analysis', 'Figures');

NbSubject = size(SubjectList,1)

NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

ROI_Coverage=nan(NbSubject,size(ListSubROI,1));
ROI_Size=nan(size(ListSubROI,1),NbLayers+1,NbSubject)

Fontsize = 8;

Jitter = linspace(0.015,0.55,NbSubject);
Jitter = Jitter(randperm(NbSubject));

%% Get Data
for SubjInd = 1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    AnalysisFolder = fullfile(SubjectFolder, 'Analysis', 'ROI');
    
    for ROI_Ind=1:size(ROIs,1)
        
        cd(fullfile(AnalysisFolder, ROIs{ROI_Ind,1}, TargetLayerFile))
        
        load(strcat('Data_', TargetLayerFile, '.mat'), 'Data_ROI', 'VoxelCount')
        
        for iSubROI=1:size(ListSubROI,1)
            
            for i=1:length(Data_ROI)
                
                temp = find(strcmp(Data_ROI(i).SubROI,ListSubROI));
                
                if temp
                    ROI_Coverage(SubjInd,temp) = Data_ROI(i).ROI_Coverage(1)/Data_ROI(i).ROI_Coverage(2);
                    ROI_Size(temp,1:NbLayers+1,SubjInd) = VoxelCount(i,1:end);
                end
                
            end
            
        end
        
    end
    
end


%% Plot ROI coverage
cd(FigureFolder)

close all

ListSubROI=strrep(ListSubROI,'_L_MNI','');
ListSubROI=strrep(ListSubROI,'_R_MNI','');
ListSubROI=strrep(ListSubROI,'_L','');
ListSubROI=strrep(ListSubROI,'_R','');

figure('name', 'ROI coverage left', 'position', [100 100 800 1000])

hold on

for iROI = 1:5
    
    for iSubj=1:size(ROI_Coverage,1)
        plot(iROI+Jitter(iSubj), ROI_Coverage(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-.1, nanmean(ROI_Coverage(:,iROI)), nansem(ROI_Coverage(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

for iROI = 6:12
    
    for iSubj=1:size(ROI_Coverage,1)
        plot(iROI+1+Jitter(iSubj), ROI_Coverage(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-.1, nanmean(ROI_Coverage(:,iROI)), nansem(ROI_Coverage(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

tmp = char(ListSubROI');
tmp = [tmp(1:5,:) ; repmat(' ',1,size(tmp,2)); tmp(6:12,:)];

axis([0.5 size(tmp,1)+.5 0 1.1])

set(gca,'tickdir', 'out', ...
    'xtick', 1:(size(tmp,1)), ...
    'xticklabel', tmp)

print(gcf, strcat('ROI coverage Left.tif'), '-dtiff')

% right
figure('name', 'ROI coverage right', 'position', [100 100 800 1000])

hold on

for iROI = 13:17
    
    for iSubj=1:size(ROI_Coverage,1)
        plot(iROI-12+Jitter(iSubj), ROI_Coverage(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-12-.1, nanmean(ROI_Coverage(:,iROI)), nansem(ROI_Coverage(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

for iROI = 18:size(ListSubROI,1)
    
    for iSubj=1:size(ROI_Coverage,1)
        plot(iROI-11+Jitter(iSubj), ROI_Coverage(iSubj,iROI), 'o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:))
    end
    
    errorbar(iROI-11-.1, nanmean(ROI_Coverage(:,iROI)), nansem(ROI_Coverage(:,iROI)),...
        ' o', 'color', 'b', 'MarkerFaceColor', 'b')
    
end

tmp = char(ListSubROI');
tmp = [tmp(1:5,:) ; repmat(' ',1,size(tmp,2)); tmp(6:12,:)];

axis([0.5 size(tmp,1)+.5 0 1.1])

set(gca,'tickdir', 'out', ...
    'xtick', 1:(size(tmp,1)), ...
    'xticklabel', tmp)

print(gcf, strcat('ROI coverage Right.tif'), '-dtiff')


%% Plot number of voxel per layer
close all

figure('name', 'ROI size Auditory', 'position', [100 100 1400 1200])

for iROI = 1:5
    
    MAX=ROI_Size([iROI iROI+12],2:end,:);
    MAX=max(MAX(:));
    
    subplot(5,2,2*iROI-1)
    hold on
    
    for iSubj=1:size(ROI_Size,3)
        tmp=ROI_Size(iROI,2:end,iSubj);
        plot((1:NbLayers)+Jitter(iSubj), fliplr(tmp), ' o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:), 'MarkerSize', 3)
        clear tmp
    end
    
    errorbar((1:NbLayers)-.1, fliplr(nanmean(ROI_Size(iROI,2:end,:),3)), ...
        fliplr(nanstd(ROI_Size(iROI,2:end,:),3)), ' o', 'color', 'b', 'MarkerFaceColor', 'b', ...
    'MarkerSize', 3)
    
    axis([0.5 7 0 MAX+100])
    
    t=ylabel('Nb voxels');
    set(t,'fontsize',Fontsize);
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', 1:NbLayers, 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    title([strrep(ListSubROI{iROI},'_',' ') ' left'])   
    
    subplot(5,2,2*iROI)
    hold on
    
    for iSubj=1:size(ROI_Size,3)
        tmp=ROI_Size(iROI+12,2:end,iSubj);
        plot((1:NbLayers)+Jitter(iSubj), fliplr(tmp), ' o', 'color', COLOR_Subject(iSubj,:),...
            'MarkerFaceColor', COLOR_Subject(iSubj,:), 'MarkerSize', 3)
    end
    
    errorbar((1:NbLayers)-.1, fliplr(nanmean(ROI_Size(iROI+12,2:end,:),3)), ...
        fliplr(nanstd(ROI_Size(iROI+12,2:end,:),3)), ' o', 'color', 'b', 'MarkerFaceColor', 'b', ...
    'MarkerSize', 3)
    
    axis([0.5 7 0 MAX+100])
    
%     t=ylabel('Number of voxels');
%     set(t,'fontsize',Fontsize);
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', 1:NbLayers, 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    title([strrep(ListSubROI{iROI},'_',' ') ' right'])
    
end

subplot(5,2,2*iROI-1)
t=xlabel('Layers');
set(t,'fontsize',Fontsize);

subplot(5,2,2*iROI)
t=xlabel('Layers');
set(t,'fontsize',Fontsize);

print(gcf, strcat('ROI size Auditory.tif'), '-dtiff')


figure('name', 'ROI size Visual', 'position', [100 100 1200 1000])

for iROI = 6:12
    
    MAX=ROI_Size([iROI iROI+12],2:end,:);
    MAX=max(MAX(:));
    
    subplot(7,2,2*iROI-11)
    hold on
    
    for iSubj=1:size(ROI_Size,3)
        tmp=ROI_Size(iROI,2:end,iSubj);
        plot((1:NbLayers)+Jitter(iSubj), fliplr(tmp), ' o', 'color', COLOR_Subject(iSubj,:), ...
            'MarkerFaceColor', COLOR_Subject(iSubj,:), 'MarkerSize', 3)
        clear tmp
    end
    
    errorbar((1:NbLayers)-.1, fliplr(nanmean(ROI_Size(iROI,2:end,:),3)), ...
        fliplr(nanstd(ROI_Size(iROI,2:end,:),3)), ' o', 'color', 'b', 'MarkerFaceColor', 'b', ...
    'MarkerSize', 3)
    
    axis([0.5 7 0 MAX+100])
    
    t=ylabel('Nb voxels');
    set(t,'fontsize',Fontsize);
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', 1:NbLayers, 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    title([strrep(ListSubROI{iROI},'_',' ') ' left'])
    
    
    
    subplot(7,2,2*iROI-10)
    hold on
    
    for iSubj=1:size(ROI_Size,3)
        tmp=ROI_Size(iROI+12,2:end,iSubj);
        plot((1:NbLayers)+Jitter(iSubj), fliplr(tmp), ' o', 'color', COLOR_Subject(iSubj,:),...
            'MarkerFaceColor', COLOR_Subject(iSubj,:), 'MarkerSize', 3)
    end
    
    errorbar((1:NbLayers)-.1, fliplr(nanmean(ROI_Size(iROI+12,2:end,:),3)), ...
        fliplr(nanstd(ROI_Size(iROI+12,2:end,:),3)), ' o', 'color', 'b', 'MarkerFaceColor', 'b', ...
    'MarkerSize', 3)
    
    axis([0.5 7 0 MAX+100])
    
%     t=ylabel('Number of voxels');
%     set(t,'fontsize',Fontsize);
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', 1:NbLayers, 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    title([strrep(ListSubROI{iROI},'_',' ') ' right'])
    
end

subplot(7,2,2*iROI-11)
t=xlabel('Layers');
set(t,'fontsize',Fontsize);

subplot(7,2,2*iROI-10)
t=xlabel('Layers');
set(t,'fontsize',Fontsize);

print(gcf, strcat('ROI size Visual.tif'), '-dtiff')

%%
cd(StartFolder)