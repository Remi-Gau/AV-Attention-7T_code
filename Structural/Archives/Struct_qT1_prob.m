%% Get T1 intensities in the ROI, TE 1.0, TE 1.1 and their intersections
clear; clc

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)

NbLayers = 10;

ProbaThres = .1;

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

MaskOri.ROI(1) = struct('name', 'TE1.0', 'fname', 'rrwmAuditory_Te10.nii');
MaskOri.ROI(end+1) = struct('name', 'TE1.1', 'fname', 'rrwmAuditory_Te11.nii');
MaskOri.ROI(end+1) = struct('name', 'TE1.2', 'fname', 'rrwmAuditory_Te12.nii');


for SubjInd = 1:size(SubjectList,1)
    
    Mask = MaskOri;
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    PMap_Folder = fullfile(SubjectFolder, 'PMap');
    
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    LayerFile = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', 'T1_10_Layers.nii');
    
    T1File = fullfile(SubjectFolder, 'Structural', 'CBS', ...
        ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii']);
    

    LayersVol = spm_read_vols(spm_vol(LayerFile));
    
    T1Vol = spm_read_vols(spm_vol(T1File));
    
    cd(PMap_Folder)
    
    try
    ROI_VOL = spm_read_vols(spm_vol(char({Mask.ROI.fname}')))>ProbaThres;
    catch
        for iImg=1:numel(Mask.ROI) 
            gunzip([Mask.ROI(iImg).fname '.gz'])
        end
        ROI_VOL = spm_read_vols(spm_vol(char({Mask.ROI.fname}')))>ProbaThres;
    end
    
    ROI_VOL(:,:,:,end+1) = any(ROI_VOL,4);
    Mask.ROI(end+1) = struct('name', 'TE', 'fname', '');
    
    qT1_Profile = nan(NbLayers,1);
    
    for LayerInd=1:NbLayers
        
        VoxPerLayer(LayerInd) = numel(find(LayersVol==LayerInd));
        
    end
    
    disp(VoxPerLayer)

    for ROI_Ind = 1:size(ROI_VOL,4)
        
        tmp = ROI_VOL(:,:,:,ROI_Ind);

        for LayerInd=1:NbLayers
            
            tmp2 = intersect(find(LayersVol==LayerInd), find(tmp));
            
            VoxPerLayer(LayerInd) = numel(tmp2);
            
            qT1_Profile(LayerInd,1) = mean(T1Vol(tmp2));
            qT1_Profile(LayerInd,2) = std(T1Vol(tmp2));
            
            clear tmp2
        end
        
        save(fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', Mask.ROI(ROI_Ind).name, ...
            '_Prob', num2str(ProbaThres), ...
            '_', num2str(NbLayers), '.mat')), 'qT1_Profile', 'VoxPerLayer')
        
        clear tmp
    end
    clear ROIs qT1_Profile
    
    
end
clear SubjInd



