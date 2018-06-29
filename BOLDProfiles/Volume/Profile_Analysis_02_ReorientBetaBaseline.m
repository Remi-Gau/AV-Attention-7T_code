clear; clc;

% StartDirectory = fullfile(pwd, '..','..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';

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
    
    matlabbatch = {};
    Files2Reorient = {};
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    TransferFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Transfer');
    
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_baseline', 'Beta');
    
    
    %%
    cd(AnalysisFolder)

    BetaFiles = dir('beta*.nii');
    for i=1:numel(BetaFiles)
        Files2Reorient{i,1} = fullfile(AnalysisFolder, [BetaFiles(i).name ',1']); %#ok<*SAGROW>
    end
    
    Files2Reorient{end+1,1} = fullfile(AnalysisFolder, 'mask.nii,1');
    
    MeanFile = dir('meanUR*.nii');
    Files2Reorient{end+1,1} = fullfile(AnalysisFolder, [MeanFile.name ',1']);

    
    %%
    cd(TransferFolder)
    
    ReorientFiles = dir('ReorientMatrix_*.mat');
    for iFile = 1:numel(ReorientFiles)
        load(fullfile(TransferFolder, ReorientFiles(iFile).name))
        matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient;
        matlabbatch{end}.spm.util.reorient.transform.transM = M;
        matlabbatch{end}.spm.util.reorient.prefix = '';
    end
    
    CoregFiles = dir('CoregMatrix_*.mat');
    for iFile = 1:numel(CoregFiles)
        load(fullfile(TransferFolder, CoregFiles(iFile).name))
        matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient;
        matlabbatch{end}.spm.util.reorient.transform.transM = M;
        matlabbatch{end}.spm.util.reorient.prefix = '';
    end

    %%
    cd(AnalysisFolder)
    save (strcat('ReorientBeta_Subject_', SubjID, '_jobs.mat'));
    
    spm_jobman('run', matlabbatch)

end