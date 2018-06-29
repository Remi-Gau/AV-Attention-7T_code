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
    
    matlabbatch = {};
    
    SubjID = SubjectList(SubjInd,:);
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'FFX_Block_Smooth');
    TransferFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Transfer');
    
    cd(AnalysisFolder)
    Files2Reorient = {};
    MaskFiles = dir('spmT*.nii');
    for i=1:numel(MaskFiles)
        Files2Reorient{i,1} = fullfile(AnalysisFolder, [MaskFiles(i).name ',1']); %#ok<SAGROW>
    end
    
    disp(Files2Reorient)
    
    cd(TransferFolder)
    
    ReorientFiles = dir('ReorientMatrix_*.mat');
    
    ReorientFiles.name
    
    for iFile = 1:numel(ReorientFiles)
        load(fullfile(TransferFolder, ReorientFiles(iFile).name))
        matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient; %#ok<*SAGROW>
        matlabbatch{end}.spm.util.reorient.transform.transM = M;
        matlabbatch{end}.spm.util.reorient.prefix = '';
    end

    
    CoregFiles = dir('CoregMatrix_*.mat');
    
    CoregFiles.name
    
    for iFile = 1:numel(CoregFiles)
        load(fullfile(TransferFolder, CoregFiles(iFile).name))
        matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient;
        matlabbatch{end}.spm.util.reorient.transform.transM = M;
        matlabbatch{end}.spm.util.reorient.prefix = '';
    end
    
    save (fullfile(TransferFolder,strcat('ReorientActMask_Subject_', SubjID, '_jobs.mat')));
    
    spm_jobman('run', matlabbatch)
    
end