%%
clc; clear;

MinLayer = 1;
NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];


StartFolder=fullfile(pwd, '..', '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

Mask.ROI(1) = struct('name', 'A1_V_Act_surf', 'fname', 'A1_V_act.nii');
Mask.ROI(end+1) = struct('name', 'PT_V_Act_surf', 'fname', 'PT_V_act.nii');
Mask.ROI(end+1) = struct('name', 'V1_A_Act_surf', 'fname', 'V1_A_act.nii');
Mask.ROI(end+1) = struct('name', 'V23_A_Act_surf', 'fname', 'V2-3_A_act.nii');

Mask.ROI(end+1) = struct('name', 'A1_V_Deact_surf', 'fname', 'A1_V_deact.nii');
Mask.ROI(end+1) = struct('name', 'PT_V_Deact_surf', 'fname', 'PT_V_deact.nii');
Mask.ROI(end+1) = struct('name', 'V1_A_Deact_surf', 'fname', 'V1_A_deact.nii');
Mask.ROI(end+1) = struct('name', 'V23_A_Deact_surf', 'fname', 'V2-3_A_deact.nii');

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention',...
    'Baseline'};


for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2', 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'FFX_Baseline');
    

    Data_Folder = fullfile(SubjectFolder, 'Transfer');
    ROI_Folder = fullfile(Data_Folder, 'ROI');
    
    BackUpFolder = fullfile(SubjectFolder, 'FFX_baseline', 'Beta');
    
    
    
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Volumes', TargetLayerFile);
    [~,~,~] = mkdir(Results_Folder);
    
    
    
    %% Gets ROIs
    
    % Gets global mask from GLM and ROI masks for the data
    Mask.global.hdr = spm_vol(fullfile(BackUpFolder, 'rmask.nii'));
    Mask.global.img = logical(spm_read_vols(Mask.global.hdr));


    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            Mask.ROI(iROI).hdr = spm_vol(fullfile(ROI_Folder, Mask.ROI(iROI).fname));
        end
    end
    
    hdr = cat(1, Mask.ROI.hdr);
    sts = spm_check_orientations([Mask.global.hdr; hdr]);
    if sts ~= 1
        error('Images not in same space!');
    end
    clear sts hdr i
    
    % Create mask in XYZ format (both world and voxel coordinates)
    [X, Y, Z] = ind2sub(size(Mask.global.img), find(Mask.global.img));
    Mask.global.XYZ = [X'; Y'; Z']; % XYZ format
    clear X Y Z
    Mask.global.size = size(Mask.global.XYZ, 2);
    Mask.global.XYZmm = Mask.global.hdr.mat(1:3,:) ...
        * [Mask.global.XYZ; ones(1, Mask.global.size)]; % voxel to world transformation
    
    % Combine masks
    fprintf(' Open ROI images\n')
    xY.def = 'mask';
    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            xY.spec = fullfile(ROI_Folder, Mask.ROI(iROI).fname);
            [xY, Mask.ROI(iROI).XYZmm, j] = spm_ROI(xY, Mask.global.XYZmm);
            Mask.ROI(iROI).XYZ = Mask.global.XYZ(:,j);
            Mask.ROI(iROI).size = size(Mask.ROI(iROI).XYZ, 2);
            A = spm_read_vols(Mask.ROI(iROI).hdr);
            A(isnan(A)) = 0;
            Mask.ROI(iROI).size(2) = sum(logical(A(:)));
            clear A
        end
    end
    clear xY j i
    
    
    %% Open layer file and gets indices of voxels in each layer
    
    fprintf(' Open layer label images & identify layer label for each voxel\n')
    
    LayerVol = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', ...
        [TargetLayerFile '.nii']);
    
    LayerVolHdr = spm_vol(LayerVol);
    
    % Number of voxel of each ROI and intersection of ROI
    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            Mask.ROI(iROI).Layers = spm_get_data(LayerVolHdr, Mask.ROI(iROI).XYZ);
            for iLayer = 1:max(Mask.ROI(iROI).Layers)
                Mask.ROI(iROI).LayersVox(1,iLayer) = sum(Mask.ROI(iROI).Layers==iLayer);
            end
        end
    end
    
    clear LayerVol LayerVolHdr
    
    %% Opens beta images
    fprintf(' Identifying the relevant beta images\n')
    
   
    for CondInd = 1:length(Conditions_Names)
        
        BetaImg{CondInd,1} =  spm_select('FPList', BackUpFolder,...
            ['^' strrep(Conditions_Names{CondInd}, ' ', '') '.nii$']); %#ok<SAGROW>

    end
    
    BetaImg

    
    
    %% Opens session constant images
    
    %% AVERAGING
    
    fprintf(' Averaging for ROI:\n')
    
    for iROI=1:length(Mask.ROI)
        
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            
            clear Data_ROI
            
            fprintf(['  ' Mask.ROI(iROI).name '\n'])
            
            Data_ROI.name = Mask.ROI(iROI).name;
            Data_ROI.Voxelcount = [Mask.ROI(iROI).size(1)  Mask.ROI(iROI).LayersVox];
            spe = repmat('%i ', [1,max(Mask.ROI(iROI).Layers)+1]);
            fprintf(['   Voxel per layer ' spe '\n'], Data_ROI.Voxelcount)
            
            Data_ROI.ROI_Coverage = Mask.ROI(iROI).size;
            fprintf('   Coverage %0.2f\n', Mask.ROI(iROI).size(1)/Mask.ROI(iROI).size(2))
            
            if ~any(Data_ROI.ROI_Coverage==0) && sum(Data_ROI.Voxelcount~=0)==(NbLayers+1)
                
                
                for CondInd = 1:length(Conditions_Names) % For each Condition
                    
                    clear SourceVolume SourceVolumeStruc
                    
                    SourceVolumeStruc = spm_vol(BetaImg{CondInd});
                    SourceVolume = spm_get_data(SourceVolumeStruc, Mask.ROI(iROI).XYZ);
                    
                    for LayerInd = 0:max(Mask.ROI(iROI).Layers) % Averages over voxels of a given layer
                        
                        Vox2Sel = Mask.ROI(iROI).Layers==LayerInd;
                        
                        Data_ROI.LayerMean(LayerInd+1,CondInd) = ...
                            nanmean(SourceVolume(Vox2Sel));
                        Data_ROI.LayerMedian(LayerInd+1,CondInd) = ...
                            nanmedian(SourceVolume(Vox2Sel));
                        Data_ROI.LayerSTD(LayerInd+1,CondInd) = ...
                            nanstd(SourceVolume(Vox2Sel));
                        Data_ROI.LayerSEM(LayerInd+1,CondInd)  = ...
                            nansem(SourceVolume(Vox2Sel));
                        
                    end
                    clear tmp LayerInd SourceVolume
                    
                end

            end
            
            cd(Results_Folder)

            save(strcat('Data_Baseline_', Data_ROI.name, '_', TargetLayerFile, '.mat'), 'Data_ROI')

            clear CondInd
            
        end
        
    end
    
    
    
end


cd(StartFolder)
