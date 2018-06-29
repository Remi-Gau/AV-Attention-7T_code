clear; clc;

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


StartDirectory = pwd;

addpath(fullfile(StartDirectory, 'SubFun'))


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

% ROI={
%     'wV1_MNI.nii'     'wROI_Visual_hOC1_L_MNI.nii' 'wROI_Visual_hOC1_R_MNI.nii' '' '' '' ''; ...
%     'wV2_MNI.nii'     'wROI_Visual_hOc2_L_MNI.nii' 'wROI_Visual_hOc2_R_MNI.nii' '' '' '' ''; ...
%     'wV3_MNI.nii'     'wROI_Visual_hOc3d_L_MNI.nii' 'wROI_Visual_hOc3d_R_MNI.nii' 'wROI_Visual_hOc3v_L_MNI.nii' 'wROI_Visual_hOc3v_R_MNI.nii' '' ''; ...
%     'wV4_MNI.nii'     'wROI_Visual_hOc4d_L_MNI.nii' 'wROI_Visual_hOc4d_R_MNI.nii' 'wROI_Visual_hOc4v_L_MNI.nii' 'wROI_Visual_hOc4v_R_MNI.nii' '' ''; ...
%     'wV5_MNI.nii'     'wROI_Visual_hOc5_L_MNI.nii'  'wROI_Visual_hOc5_R_MNI.nii' '' '' '' ''; ...
%     'wTE_1.0_MNI.nii' 'wROI_auditory_Te10_L_MNI.nii' 'wROI_auditory_Te10_R_MNI.nii' '' '' '' ''; ...
%     'wTE_1.1_MNI.nii' 'wROI_auditory_Te11_L_MNI.nii' 'wROI_auditory_Te11_R_MNI.nii' '' '' '' ''; ...
%     'wTE_1.2_MNI.nii' 'wROI_auditory_Te12_L_MNI.nii' 'wROI_auditory_Te12_R_MNI.nii' '' '' '' ''; ...
%     'wTE_MNI.nii'     'wROI_auditory_Te10_L_MNI.nii' 'wROI_auditory_Te10_R_MNI.nii' 'wROI_auditory_Te11_L_MNI.nii' 'wROI_auditory_Te11_R_MNI.nii' 'wROI_auditory_Te12_L_MNI.nii' 'wROI_auditory_Te12_R_MNI.nii'};


ROI={
    'wTE_1.0_p>.6.nii'     'wTE_1.0_p>.6_lh.nii' 'wTE_1.0_p>.6_rh.nii'; ...
    'wTE_1.1_p>.6.nii'     'wTE_1.1_p>.6_lh.nii' 'wTE_1.1_p>.6_rh.nii'; ...
    'wTE_1.2_p>.6.nii'     'wTE_1.2_p>.6_lh.nii' 'wTE_1.2_p>.6_rh.nii'};

flags = struct(...
    'mask', 0, ...
    'mean', 0, ...
    'interp', 0, ...
    'which', 1, ...
    'wrap', [0 0 0], ...
    'prefix', 'rr' ...
    );


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
                            
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);                    
    
    ROI_Folder = fullfile(SubjectFolder, 'ROI_MNI');
    
    AnalysisFolder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    
%     copyfile(fullfile(ROI_Folder,'w*S1*.*'), fullfile(AnalysisFolder));
%     copyfile(fullfile(ROI_Folder,'wProbRet*d.nii'), fullfile(AnalysisFolder));
%     copyfile(fullfile(ROI_Folder,'wProbRet*v.nii'), fullfile(AnalysisFolder));
    
    copyfile(fullfile(ROI_Folder,'wTE_1.*.nii'), fullfile(AnalysisFolder));
    
    
%     copyfile(fullfile(SubjectFolder, 'FFX_NoSmooth', 'mask.nii'), fullfile(AnalysisFolder));

    %% Merge ROIs
    cd(AnalysisFolder)
    for iROI=1:size(ROI,1)
        
        clear hdr vol 
        for img=2:size(ROI,2);
            if ~strcmp(ROI{iROI,img},'')
                hdr{img-1} = ROI{iROI,img}; %#ok<SAGROW>
            else
                break
            end
        end
        
        hdr = spm_vol(char(hdr'));
        vol = spm_read_vols(hdr);
        vol = any(vol,4);
        
        hdr = hdr(1);
        
        hdr.fname = ROI{iROI,1};
        
        spm_write_vol(hdr,vol);
        
    end
    
%     delete('rr*.nii')
%     delete('wROI*.nii')

    %% Reorient ROIs
    cd(AnalysisFolder)
    
    Files2Reorient = {};
    
%     filesdir = dir('w*S1*.*');
%     tmp = dir('wProbRet*d.nii');
%     tmp2 = dir('wProbRet*v.nii');
%     filesdir = [filesdir;tmp;tmp2]; %#ok<*AGROW>
    
    filesdir = dir('wTE_1.*_p>.6.nii');  

    for iROI=1:numel(filesdir)
        Files2Reorient{end+1,1} = fullfile(AnalysisFolder,filesdir(iROI).name); %#ok<*SAGROW>
    end
    clear tmp i
    
    cd(fullfile(SubjectFolder, 'Transfer'));
    
    MatFiles = dir('*mat');
    MatFiles = char(MatFiles.name);
    
    Reorient = strcmp(cellstr(MatFiles(:,1:14)),'ReorientMatrix');
    if any(Reorient)
        if sum(Reorient)>1
            error('Images have been reoriented more than once.')
            return
        else
            load(deblank(MatFiles(Reorient,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    clear Reorient M
    
    Coregister = strcmp(cellstr(MatFiles(:,1:11)),'CoregMatrix');
    if any(Coregister)
        if sum(Coregister)>1
            error('Images have been coregistered more than once.')
            return
        else
            load(deblank(MatFiles(Coregister,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    
    clear Files2Reorient M
    
    %% Reslice
    cd(AnalysisFolder)

    Files2Reslice = {};
    matlabbatch = {};
    
%     ROIFiles = dir('w*S1*.*'); %#ok<NASGU>
%     tmp = dir('wProbRet*d.nii');
%     tmp2 = dir('wProbRet*v.nii');
%     ROIFiles = [filesdir;tmp;tmp2]; %#ok<*AGROW>
    
    ROIFiles = dir('wTE_1.*_p>.6.nii'); 
    
    for i=1:numel(ROIFiles)
        Files2Reslice{end+1,1} = fullfile(AnalysisFolder, [ROIFiles(i).name ',1']);
    end
    Files2Reslice{end+1,1} = fullfile(AnalysisFolder, 'mask.nii,1');

    matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'rr';

    save (strcat('ResliceROIs_Subject_', SubjID, '_jobs.mat'));
    
    spm_jobman('run', matlabbatch)


end