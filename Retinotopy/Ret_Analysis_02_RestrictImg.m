clc; clear;

StartDirectory = pwd;

addpath(fullfile(StartDirectory, 'SubFun'))

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

 for SubjInd = 1:size(SubjectList,1)
     
     SubjID = SubjectList(SubjInd,:);
     
     SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], ...
         'Retinotopy', 'MIPAV');
     
     fprintf('\nAnalysing subject %s\n', SubjID)
         
     AnalysisFolder = fullfile(SubjectFolder);
     cd(AnalysisFolder)
     BetaFiles = dir('rPolar*.nii');
     for iBeta = 1:numel(BetaFiles)
         hdr = spm_vol(BetaFiles(iBeta).name);
         vol = spm_read_vols(hdr);
         vol(isnan(vol))=0;
         spm_write_vol(hdr,vol);
         hdr=[]; vol=[];
     end
     
 end
