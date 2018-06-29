%% Get T1 intensities in the ROI, TE 1.0, TE 1.1 and their intersections
clear; clc

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)

NbLayers = 10;

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

for SubjInd = 1:size(SubjectList,1)
    
    clear Mask
    
%     Mask.ROI(1) = struct('name', 'V1', 'fname', 'V1_surf.nii');
%     Mask.ROI(end+1) = struct('name', 'V2-3_surf', 'fname', 'V2-3_surf.nii');
%     Mask.ROI(end+1) = struct('name', 'A1_surf', 'fname', 'A1_surf.nii');
%     Mask.ROI(end+1) = struct('name', 'PT_BT', 'fname', 'PT_surf.nii');
    
    
%     Mask.ROI(1) = struct('name', 'V1_surf_thres', 'fname', 'V1_surf_thres.nii');
    Mask.ROI(1) = struct('name', 'V2-3_surf_thres', 'fname', 'V2-3_surf_thres.nii');
%     Mask.ROI(1) = struct('name', 'PT_surf_thres', 'fname', 'PT_surf_thres.nii');
%     Mask.ROI(1) = struct('name', 'A1_surf', 'fname', 'A1_surf.nii');

    
%     Mask.ROI(end+1) = struct('name', 'TE1.0_surf', 'fname', 'TE10_surf.nii');
%     Mask.ROI(end+1) = struct('name', 'TE1.1_surf', 'fname', 'TE11_surf.nii');
%     Mask.ROI(end+1) = struct('name', 'TE1.2_surf', 'fname', 'TE12_surf.nii');
    
    
    % Mask.ROI(1) = struct('name', 'TE', 'fname', 'rrwTE_MNI.nii');
    
    % Mask.ROI(end+1) = struct('name', 'pSTG', 'fname', 'rrwSTG_Post_AAL.nii');
    
%     Mask.ROI(1) = struct('name', 'TE1.0', 'fname', 'rrwTE_1.0_MNI.nii');
    % Mask.ROI(end+1) = struct('name', 'TE1.1', 'fname', 'rrwTE_1.1_MNI.nii');
    % Mask.ROI(end+1) = struct('name', 'TE1.2', 'fname', 'rrwTE_1.2_MNI.nii');
    %
    % Mask.ROI(end+1) = struct('name', 'TE_surf', 'fname', 'TE_surf.nii');
    %
    % Mask.ROI(end+1) = struct('name', 'pSTG_surf', 'fname', 'pSTG_surf.nii');
    
    %
    % Mask.ROI(end+1) = struct('name', 'V1', 'fname', 'rrwProbRet_V1.nii');
    % Mask.ROI(end+1) = struct('name', 'V2', 'fname', 'rrwProbRet_V2.nii');
    % Mask.ROI(end+1) = struct('name', 'V3', 'fname', 'rrwProbRet_V3.nii');
    %
    % Mask.ROI(end+1) = struct('name', 'A1', 'fname', 'A1.nii');
    %
    % Mask.ROI(end+1) = struct('name', 'TE1.0-p>.6', 'fname', 'rrwTE_1.0_p>.6.nii');
    % Mask.ROI(end+1) = struct('name', 'TE1.1-p>.6', 'fname', 'rrwTE_1.1_p>.6.nii');
    % Mask.ROI(end+1) = struct('name', 'TE1.2-p>.6', 'fname', 'rrwTE_1.2_p>.6.nii');
    
    
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    Data_Folder = fullfile(SubjectFolder, 'Transfer');
    
    ROI_Folder = fullfile(Data_Folder, 'ROI');
    
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    LayerFile = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', 'T1_10_Layers.nii');
    
%     T1GLMFile = fullfile(SubjectFolder, 'Structural', 'CBS', ...
%         ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii']);
    
    T1GLMFile = fullfile(StartFolder, 'Subjects_Data', ...
        ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_roiglm_approx.nii']);
    gunzip([T1GLMFile '.gz'])

    LayersVol = spm_read_vols(spm_vol(LayerFile));
    
    T1Vol = spm_read_vols(spm_vol(T1GLMFile));
    
    cd(ROI_Folder)
    ROI_VOL = logical(spm_read_vols(spm_vol(char({Mask.ROI.fname}'))));
    
    qT1_Profile = nan(NbLayers,1);
    
    for LayerInd=1:NbLayers
        
        VoxPerLayer(LayerInd) = numel(find(LayersVol==LayerInd));
        
    end
    
    VoxPerLayer

    for ROI_Ind = 1:size(ROI_VOL,4)
        
        tmp = ROI_VOL(:,:,:,ROI_Ind);

        for LayerInd=1:NbLayers
            
            tmp2 = intersect(find(LayersVol==LayerInd), find(tmp));
            
            tmpLR = find(tmp);
            [X,Y,Z] = ind2sub(size(tmp),tmpLR);
            tmpL = intersect(find(LayersVol==LayerInd),tmpLR(X<size(tmp,1)/2));
            tmpR = intersect(find(LayersVol==LayerInd),tmpLR(X>size(tmp,1)/2));
            
            VoxPerLayer(1,LayerInd) = numel(tmp2);
            VoxPerLayerL(1,LayerInd) = numel(tmpL);
            VoxPerLayerR(1,LayerInd) = numel(tmpR);
            
            qT1_Profile(LayerInd,1) = mean(T1Vol(tmp2));
            qT1_Profile(LayerInd,2) = std(T1Vol(tmp2));
            
            qT1_ProfileL(LayerInd,1) = mean(T1Vol(tmpL));
            qT1_ProfileL(LayerInd,2) = std(T1Vol(tmpL));
            
            qT1_ProfileR(LayerInd,1) = mean(T1Vol(tmpR));
            qT1_ProfileR(LayerInd,2) = std(T1Vol(tmpR));
            
            
            clear tmp2
        end
        
        save(fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', strrep(Mask.ROI(ROI_Ind).name,'_surf',''), '_GLM_', num2str(NbLayers), '.mat')), ...
            'qT1_Profile', 'VoxPerLayer', 'qT1_ProfileL', 'VoxPerLayerL', 'qT1_ProfileR', 'VoxPerLayerR')
        
        clear tmp
    end
    clear ROIs qT1_Profile
    
    
end
clear SubjInd



