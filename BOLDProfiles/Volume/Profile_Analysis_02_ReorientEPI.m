clear; clc;

StartDirectory = fullfile(pwd, '..', '..', '..');

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


for SubjInd = 1 %1:size(SubjectList,1)
    
    Files2Reorient = {};
    
    SubjID = SubjectList(SubjInd,:);
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Nifti');
    TransferFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Transfer');
    
    cd(AnalysisFolder)
    
    DirLs = dir('0*');
    NbDir = numel(DirLs);
    
    for iDir=1:NbDir
        cd(fullfile(AnalysisFolder, DirLs(iDir).name) )
        EPIFile = dir('URS*.nii');
        EPIHdr = spm_vol(EPIFile.name);
        for iImg=1:numel(EPIHdr)
            Files2Reorient{end+1,1} = ...
                fullfile(AnalysisFolder, DirLs(iDir).name, [EPIHdr(iImg).fname ',' num2str(iImg)]); %#ok<SAGROW>
        end
    end
    

    cd(TransferFolder);
    M = dir('ReorientMatrix_*.mat');
    load(fullfile(TransferFolder, M.name))
    
    matlabbatch{1}.spm.util.reorient.srcfiles = Files2Reorient;
    matlabbatch{1}.spm.util.reorient.transform.transM = M;
    matlabbatch{1}.spm.util.reorient.prefix = '';
    
    M = dir('CoregMatrix_*.mat');
    load(fullfile(TransferFolder, M.name))

    matlabbatch{2}.spm.util.reorient.srcfiles = Files2Reorient;
    matlabbatch{2}.spm.util.reorient.transform.transM = M;
    matlabbatch{2}.spm.util.reorient.prefix = '';

    cd(AnalysisFolder)
    save (strcat('ReorientEPI_Subject_', SubjID, '_jobs.mat'));
    
%     spm_jobman('run', matlabbatch)

end