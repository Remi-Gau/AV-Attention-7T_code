% segment the T1w for each subject


clear; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults

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

if ~matlabpool('SIZE')
    matlabpool('OPEN', 3) ;
end

parfor SubjInd = 1:size(SubjectList,1)
    
    tic
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    StructuralScan = fullfile(SubjectFolder, 'FFX', 'Structural', 'UNI.nii');
    
    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    matlabbatch = {};
    
    matlabbatch{1,1}.spm.spatial.preproc.channel.vols = {StructuralScan};
    matlabbatch{1,1}.spm.spatial.preproc.channel.biasreg = 1.0000e-03;
    matlabbatch{1,1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1,1}.spm.spatial.preproc.channel.write = [1 1];
    
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(1).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,1'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(2).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,2'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(3).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,3'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(4).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,4'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(5).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,5'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.tissue(6).tpm = {'/home/SHARED/Program/Matlab/spm12/tpm/TPM.nii,6'};
    matlabbatch{1,1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1,1}.spm.spatial.preproc.tissue(6).native = [1 0];
    matlabbatch{1,1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    
    matlabbatch{1,1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1,1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1,1}.spm.spatial.preproc.warp.reg = [0 1.0000e-03 0.5000 0.0500 0.2000];
    matlabbatch{1,1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1,1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1,1}.spm.spatial.preproc.warp.samp = 2;
    matlabbatch{1,1}.spm.spatial.preproc.warp.write = [1 1];

    
    cd (fullfile(SubjectFolder, 'FFX', 'Structural'))
    
    save (strcat('Segment_', SubjID, '_matlabbatch'));
    
    fprintf('\n\n')
    disp('%%%%%%%%%%%%%%%')
    disp('    SEGMENT    ')
    disp('%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
    toc
    
end

matlabpool('CLOSE')