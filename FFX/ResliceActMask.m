clear; clc;


% StartDirectory = fullfile(pwd, '..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';

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
    
    SubjID = SubjectList(SubjInd,:);
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'FFX_Block_Smooth');
  
    cd(AnalysisFolder)
    
    Files2Reslice = {};
    matlabbatch = {};
    
    MaskFiles = dir('spm*0.1.nii');
    for i=1:numel(MaskFiles)
        Files2Reslice{i,1} = fullfile(AnalysisFolder, [MaskFiles(i).name ',1']); %#ok<*SAGROW>
    end
    
    Files2Reslice
    
        matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    save (fullfile(AnalysisFolder,strcat('ReSliceActMask_Subject_', SubjID, '_jobs.mat')));

    spm_jobman('run', matlabbatch)


end