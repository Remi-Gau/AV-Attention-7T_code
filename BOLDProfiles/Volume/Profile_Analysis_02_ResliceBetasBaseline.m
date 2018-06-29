clear; clc;

% StartDirectory = fullfile(pwd, '..','..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';
cd (StartDirectory)

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


for SubjInd = 2:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_baseline', 'Beta');
    
    cd(AnalysisFolder)
    Files2Reslice = {};
    matlabbatch = {};
    
    MeanFile = dir('meanUR*.nii');
    for i=1:numel(MeanFile)
        Files2Reslice{i,1} = fullfile(AnalysisFolder, [MeanFile(i).name ',1']); %#ok<*SAGROW>
    end
    
    BetaFiles = dir('beta*.nii');
    for i=1:numel(BetaFiles)
        Files2Reslice{i,1} = fullfile(AnalysisFolder, [BetaFiles(i).name ',1']); %#ok<*SAGROW>
    end
    
    matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    Files2Reslice
    
    Files2Reslice = {};
    Files2Reslice{1,1} = fullfile(AnalysisFolder, 'mask.nii,1');
    
    matlabbatch{2}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    matlabbatch{2}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    Files2Reslice
    
    save (strcat('ResliceMean_Subject_', SubjID, '_jobs.mat'));
    
    tic
    spm_jobman('run', matlabbatch)
    toc
    
    
end