clear; clc;

% StartDirectory = fullfile(pwd, '..','..', '..');
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
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID] ...
        , 'Transfer2');
    
    cd(AnalysisFolder)
    
    Files2Reorient = {};
    
%     BetaFiles = dir('beta*.nii');
%     for i=1:numel(BetaFiles)
%         Files2Reorient{i,1} = fullfile(AnalysisFolder, [BetaFiles(i).name ',1']); %#ok<*SAGROW>
%     end
    
%     ConFiles = dir('con*.nii');
%     for i=1:numel(ConFiles)
%         Files2Reorient{end+1,1} = fullfile(AnalysisFolder, [ConFiles(i).name ',1']);
%     end

    MeanFile = dir('mean*.nii');
    Files2Reorient{1,1} = fullfile(AnalysisFolder, [MeanFile(1).name ',1']);
    Files2Reorient{2,1} = fullfile(AnalysisFolder, [MeanFile(2).name ',1']);
    
%     Files2Reorient{end+1,1} = fullfile(AnalysisFolder, 'mask.nii,1');

    ReorientFiles = dir('ReorientMatrix_*.mat');
    for iFile = 1:numel(ReorientFiles)
        load(fullfile(AnalysisFolder, ReorientFiles(iFile).name))
        matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient;
        matlabbatch{end}.spm.util.reorient.transform.transM = M;
        matlabbatch{end}.spm.util.reorient.prefix = '';
    end
    
%     CoregFiles = dir('CoregMatrix_*.mat');
%     for iFile = 1:numel(CoregFiles)
%         load(fullfile(AnalysisFolder, CoregFiles(iFile).name))
%         matlabbatch{end+1}.spm.util.reorient.srcfiles = Files2Reorient;
%         matlabbatch{end}.spm.util.reorient.transform.transM = M;
%         matlabbatch{end}.spm.util.reorient.prefix = '';
%     end

%     save (strcat('ReorientBeta_Subject_', SubjID, '_jobs.mat'));

    save (strcat('ReorientMean_Subject_', SubjID, '_jobs.mat'));
    
    spm_jobman('run', matlabbatch)

end