clear all; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

%  Folders definitions
RootFolder = pwd;

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


for SubjInd = 2:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    cd(fullfile(SubjectFolder, 'Nifti', 'NoMoCo', '01'));
    MeanImage = dir('mean*.nii');
    MeanImage = fullfile(pwd, strcat(MeanImage.name, ',1'));
    
    ROI_Folder = fullfile(SubjectFolder, 'ROI_MNI');
    cd(ROI_Folder)
    filesdir = dir('w*.nii');
    for ImageInd=1:length(filesdir)
        ImagesFiles2Process{ImageInd,1} = fullfile(pwd, strcat(filesdir(ImageInd).name, ',1'));
    end

    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    matlabbatch{1}.spm.spatial.coreg.write.ref = {MeanImage};
    matlabbatch{1}.spm.spatial.coreg.write.source = ImagesFiles2Process;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end