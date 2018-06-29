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
    
    SubjID = SubjectList(SubjInd,:);
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Nifti');
    TransferFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Transfer');
    
    cd(AnalysisFolder)
    
    Files2Reslice = {};
    matlabbatch = {};
    
    DirLs = dir('0*');
    NbDir = numel(DirLs);
    
    for iDir=1:NbDir
        cd(fullfile(AnalysisFolder, DirLs(iDir).name) )
        EPIFile = dir('URS*.nii');
        EPIHdr = spm_vol(EPIFile.name);
        for iImg=1:numel(EPIHdr)
            Files2Reslice{end+1,1} = ...
                fullfile(AnalysisFolder, DirLs(iDir).name, [EPIHdr(iImg).fname ',' num2str(iImg)]); %#ok<SAGROW>
        end
    end 

    matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    cd(AnalysisFolder)
    save (strcat('ResliceEPI_Subject_', SubjID, '_jobs.mat'));

    spm_jobman('run', matlabbatch)


end