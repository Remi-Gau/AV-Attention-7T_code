%%
clear all; clc

StartDirectory = pwd;

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

ROIs = {...
'wROI_auditory_Te10_L_MNI.nii', 1 ; ...
'wROI_auditory_Te11_L_MNI.nii', 2 ; ...
'wROI_auditory_Te12_L_MNI.nii', 3 ; ...
'wROI_auditory_TI1_L_MNI.nii', 4 ; ...

'wROI_Visual_hOC1_L_MNI.nii', 11 ; ...
'wROI_Visual_hOc2_L_MNI.nii', 12 ; ...
'wROI_Visual_hOc3d_L_MNI.nii', 13 ; ...
'wROI_Visual_hOc3v_L_MNI.nii', 13.5 ; ...
'wROI_Visual_hOc4d_L_MNI.nii', 14 ; ...
'wROI_Visual_hOc4v_L_MNI.nii', 14.5 ; ...
'wROI_Visual_hOc5_L_MNI.nii', 15 ; ...

'wROI_auditory_Te10_R_MNI.nii', 101 ; ...
'wROI_auditory_Te11_R_MNI.nii', 102 ; ...
'wROI_auditory_Te12_R_MNI.nii', 103 ; ...
'wROI_auditory_TI1_R_MNI.nii', 104 ; ...

'wROI_Visual_hOC1_R_MNI.nii', 111 ; ...
'wROI_Visual_hOc2_R_MNI.nii', 112 ; ...
'wROI_Visual_hOc3d_R_MNI.nii', 113 ; ...
'wROI_Visual_hOc3v_R_MNI.nii', 113.5 ; ...
'wROI_Visual_hOc4d_R_MNI.nii', 114 ; ...
'wROI_Visual_hOc4v_R_MNI.nii', 114.5 ; ...
'wROI_Visual_hOc5_R_MNI.nii', 115 ; ...
};

flags = struct(...
    'sep', [2 1], ...
    'params',  [0 0 0  0 0 0], ...
    'cost_fun', 'nmi', ...
    'tol', [repmat(0.001, 1, 3), repmat(0.0005, 1, 3), repmat(0.005, 1, 3), repmat(0.0005, 1, 3)], ...
    'fwhm', [7,7], ...
    'graphics', ~spm('CmdLine'));

%  Root folder definition
StartDirectory = pwd;

cd(StartDirectory)

for SubjInd = 1:size(SubjectList,1)
    %% Subject's Identity and folders
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    AnalysisFolder = fullfile(SubjectFolder, 'Structural', 'CBS', 'ROI_MNI');
    
    
    %% Move files   
    copyfile(fullfile(SubjectFolder, 'FFX', 'Structural', 'UNI.nii'), AnalysisFolder);
    copyfile(fullfile(SubjectFolder, 'ROI_MNI', 'wROI*L_MNI.nii'), AnalysisFolder);
    copyfile(fullfile(SubjectFolder, 'ROI_MNI', 'wROI*R_MNI.nii'), AnalysisFolder);
    

    %% Creates one volume summing up all the  ROIs
    cd(AnalysisFolder);
    
    delete CoregMat*.mat

    tmp = spm_vol(ROIs{1,1});
    HDR = struct(...
        'fname',   'ALL_ROIs.nii',...
        'dim',     tmp.dim,...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     tmp.mat,...
        'pinfo',   [1 0 0]',...
        'descrip', 'All ROIs');
    
    VOL = zeros(tmp.dim);
    clear tmp
    
    for iROI = 1:size(ROIs,1)
        
        VolROI = logical(spm_read_vols(spm_vol(ROIs{iROI,1})));
        
        if any(VolROI(VolROI)==logical(VOL(VolROI)))
            tmp = VOL(VolROI);
            tmp = tmp(tmp~=0);
            warning('ROI overlap: %i voxels', length(tmp))

        end
       
        VOL(VolROI) = VOL(VolROI)+ROIs{iROI,2};
        
        clear tmp
        
    end
    
    spm_write_vol(HDR, VOL);
    
    delete w*.nii
    
    
    %% Moves UNI used for segmentation to roughly the same space as the T1 structural
    T1_Hdr = spm_vol(fullfile(AnalysisFolder,'UNI.nii'));
    T1_Vol = spm_read_vols(T1_Hdr);
    
    P_translation = zeros(12,1);
    
    P_translation(1) = 30;
    P_translation(2) = 15;
    P_translation(3) = -190;
    
    TranslationMat = spm_matrix(P_translation, 'T');
    
    T1_Hdr.mat = TranslationMat * T1_Hdr.mat;
    
    spm_write_vol(T1_Hdr, T1_Vol);
    
    
    %% Ditto ROI image
    Hdr = spm_vol(fullfile(AnalysisFolder,'ALL_ROIs.nii'));
    Vol = spm_read_vols(Hdr);
    
    P_translation = zeros(12,1);
    
    Hdr.mat = TranslationMat * Hdr.mat;
    
    spm_write_vol(Hdr, Vol);


    %% Coregister
    TargetScan = fullfile(SubjectFolder, 'Structural', 'CBS', ...
        ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii']); 
    SourceScan = fullfile(AnalysisFolder,'UNI.nii');
    
    spm_coreg_reorient_save(TargetScan,SourceScan,{fullfile(AnalysisFolder,'ALL_ROIs.nii')},flags);
    
    
    %% Reslice
    matlabbatch = {};

    matlabbatch{1}.spm.spatial.coreg.write.ref{1} = TargetScan;
    matlabbatch{1}.spm.spatial.coreg.write.source = {fullfile(AnalysisFolder,'ALL_ROIs.nii')};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

    save(fullfile(AnalysisFolder, 'ResliceBatch.mat'), 'matlabbatch')

    spm_jobman('run',matlabbatch)

    
    cd (StartDirectory)
end