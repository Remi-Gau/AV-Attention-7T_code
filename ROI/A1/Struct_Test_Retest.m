clear; close all; clc;

StartFolder = pwd;

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';...
    %     '07';... Problm with vtk files
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '14';...
    %     '15';...
    '16'
    ];

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    fprintf(['\n\nProcessing subject : ' SubjID '\n'])
    
    
    % Volumes
    ROI_L = logical(spm_read_vols(spm_vol(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1_L.nii'))));
    ROI_L(:,:,:,2) = logical(spm_read_vols(spm_vol(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1_2_L.nii'))));
    
    ROI_R = logical(spm_read_vols(spm_vol(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1_R.nii'))));
    ROI_R(:,:,:,2) = logical(spm_read_vols(spm_vol(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1_2_R.nii'))));
    
    
    CorticalMask = logical(spm_read_vols(spm_vol(fullfile(SubjectFolder, 'Structural', ...
        'CBS', 'Layering', 'T1_04_Layers.nii'))));
    
    
    % Surfaces
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_Mapping'))
    LogFileList = dir(strcat(['A1_T1_' SubjID '*lcr*inf*.vtk']));
    [~, ~, ROI_mapping_L(:,1)] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    
    LogFileList = dir(strcat(['A1_2_T1_' SubjID '*lcr*inf*.vtk']));
    [~, ~, ROI_mapping_L(:,2)] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    
    LogFileList = dir(strcat(['A1_T1_' SubjID '*rcr*inf*.vtk']));
    [~, ~, ROI_mapping_R(:,1)] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    
    LogFileList = dir(strcat(['A1_2_T1_' SubjID '*rcr*inf*.vtk']));
    [~, ~, ROI_mapping_R(:,2)] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    
    ROI_mapping_L = ROI_mapping_L==35;
    ROI_mapping_R = ROI_mapping_R==135;
    
    % VOXELS
    % LEFT
    % First delineation
    tmp = ROI_L(:,:,:,1);
    NbVox{1,1}(SubjInd,1) = sum(tmp(:));
    tmp(~CorticalMask)=0; % Restrict to cortical voxels
    NbVox{2,1}(SubjInd,1) = sum(tmp(:));
    
    % Second delineation
    tmp = ROI_L(:,:,:,2);
    NbVox{1,1}(SubjInd,2) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,1}(SubjInd,2) = sum(tmp(:));
    
    % Intersection
    tmp = all(ROI_L,4);
    NbVox{1,1}(SubjInd,3) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,1}(SubjInd,3) = sum(tmp(:));
    
    % Union
    tmp = any(ROI_L,4);
    NbVox{1,1}(SubjInd,4) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,1}(SubjInd,4) = sum(tmp(:));
    
    % RIGHT
    tmp = ROI_R(:,:,:,1);
    NbVox{1,2}(SubjInd,1) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,2}(SubjInd,1) = sum(tmp(:));
    
    tmp = ROI_R(:,:,:,2);
    NbVox{1,2}(SubjInd,2) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,2}(SubjInd,2) = sum(tmp(:));
    
    tmp = all(ROI_R,4);
    NbVox{1,2}(SubjInd,3) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,2}(SubjInd,3) = sum(tmp(:));
    
    tmp = any(ROI_R,4);
    NbVox{1,2}(SubjInd,4) = sum(tmp(:));
    tmp(~CorticalMask)=0;
    NbVox{2,2}(SubjInd,4) = sum(tmp(:));
    
    
    
    % VERTICES
    % LEFT
    % First delineation
    NbVer{1,1}(SubjInd,1) = sum(ROI_mapping_L(:,1));
    
    % Second delineation
    NbVer{1,1}(SubjInd,2) = sum(ROI_mapping_L(:,2));
    
    % Intersection
    tmp = all(ROI_mapping_L,2);
    NbVer{1,1}(SubjInd,3) = sum(tmp);
    
    % Union
    tmp = any(ROI_mapping_L,2);
    NbVer{1,1}(SubjInd,4) = sum(tmp);
    
    
    % RIGHT
    NbVer{1,2}(SubjInd,1) = sum(ROI_mapping_R(:,1));
    
    NbVer{1,2}(SubjInd,2) = sum(ROI_mapping_R(:,2));
    
    tmp = all(ROI_mapping_R,2);
    NbVer{1,2}(SubjInd,3) = sum(tmp);
    
    tmp = any(ROI_mapping_R,2);
    NbVer{1,2}(SubjInd,4) = sum(tmp);
    
    
    clear tmp ROI_L ROI_R CorticalMask ROI_mapping_R ROI_mapping_L
end

cd(StartFolder)

NbVoxBU = NbVox;

%%
close all

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

Prefix='VOX_'

for VoxVer = 1:2
    
    if VoxVer==2;
        NbVox=NbVer;
        Prefix='VER_'
    else
        NbVox=NbVoxBU;
        Prefix='VOX_'
    end;
    
    figure('name', 'T1 based A1c: test-retest', 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1]);
    
    iSubplot=1;
    
    MAX = cell2mat(NbVox);
    MAX = max(MAX(:));
    
    
    for i=1:size(NbVox,1)
        for j=1:size(NbVox,2)
            
            subplot(size(NbVox,1),size(NbVox,2),iSubplot)
            
            bar(NbVox{i,j}(:,1:3))
            
            set(gca,'tickdir', 'out', 'xtick', 1:size(SubjectList,1), 'xticklabel', SubjectList, ...
                'ticklength', [0.01 0.01], 'fontsize', Fontsize)
            
            H = axis;
            axis([0.5 size(SubjectList,1)+.5 H(3) MAX])
            
            t=ylabel('Number of voxels');
            set(t,'fontsize',Fontsize);
            
            t=xlabel('Subject');
            set(t,'fontsize',Fontsize);
            
            iSubplot=iSubplot+1;
        end
        
    end
    
    print(gcf, [Prefix 'A1_test-retest.tif'], '-dtiff')
    
    
    
    figure('name', 'T1 based A1c: test-retest correlation', 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1]);
    
    iSubplot=1;
    
    for i=1:size(NbVox,1)
        for j=1:size(NbVox,2)
            
            subplot(size(NbVox,1),size(NbVox,2),iSubplot)
            
            hold on
            
            for SubjInd=1:size(NbVox{i,j},1)
                plot(NbVox{i,j}(SubjInd,1), NbVox{i,j}(SubjInd,2), ...
                    'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
                    'MarkerFaceColor', COLOR_Subject(SubjInd,:))
            end
            [B,DEV,STATS] = glmfit(NbVox{i,j}(:,1),NbVox{i,j}(:,2),'normal');
            plot(NbVox{i,j}(:,1), NbVox{i,j}(:,1)*B(2)+B(1), 'color', [0.08 0.17 .55]);
            
            axis([0 MAX 0 MAX])
            axis square
            
            t=ylabel('Number of voxels delineation 2');
            set(t,'fontsize',Fontsize);
            
            t=xlabel('Number of voxels delineation 1');
            set(t,'fontsize',Fontsize);
            
            iSubplot=iSubplot+1;
            
        end
    end
    
    print(gcf, [Prefix 'A1_test-retest_correlation.tif'], '-dtiff')
    
    
    
    figure('name', 'INTERSECTION / UNION', 'Position', [100, 100, 1500, 1000], ...
        'Color', [1 1 1], 'Visible', 'on');
    
    hold on
    
    Absi = 1;
    
    for i=1:size(NbVox,1)
        for j=1:size(NbVox,2)
            
            tmp = NbVox{i,j}(:,3)./NbVox{i,j}(:,4);
            
            errorbar(Absi, mean(tmp), nansem(tmp), 'MarkerFaceColor',[0 0 1],...
                'MarkerEdgeColor',[0 0 1], 'Marker','o')
            for SubjInd=1:size(tmp,1)
                plot(Absi+.1+.05*SubjInd, tmp(SubjInd,1), ...
                    'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
                    'MarkerFaceColor', COLOR_Subject(SubjInd,:))
                
            end
            Absi = Absi + 1;
            
        end
    end
    
    TMP = {...
        'A1_L vs A1_2_L';...
        'A1_R vs A1_2_R';...
        'Cortical: A1_L vs A1_2_L';...
        'Cortical: A1_R vs A1_2_R'};
    
    set(gca,'tickdir', 'out', 'xtick', 1:4 , ...
        'xticklabel', char(TMP) , 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    axis([0.5 Absi 0 1])
    
    t=ylabel('INTERSECTION / UNION');
    set(t,'fontsize',Fontsize);
    
    
    print(gcf, [Prefix 'A1_test-retest_Inter_VS_Union.tif'], '-dtiff')
    
end