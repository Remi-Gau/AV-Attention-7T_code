clear all; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

%  Folders definitions
RootFolder = pwd;

%% Phase enconded lines (PELines) and ReadOutTime
% PELines = ((BaseResolution * PartialFourier)/iPat) + ((iPat-1)/iPAT) * ReferenceLines) =
% ReadoutDuration = PELines * InterEchoSpacing

% GRAPPA=iPAT4 ; Partial Fourrier=6/8 ; 60 sli ; TE=18ms ; Res=1.5 mm
% Bandwidth Per Pixel Phase Encode = 22.755


%% According to Robert Trampel

% For distortion correction: ignore Partial Fourrier and references lines
% BaseResolution/iPAT = PELines

% Effective echo spacing: 2 ways to calculate, should be the same
% 1/(Bandwidth Per Pixel Phase Encode * Reconstructed phase lines) -->   ms
% echo spacing (syngo) / iPAT

% SPM Total readout time = 1/"Bandwidth Per Pixel Phase Encode", stored in
% DICOM tag (0019, 1028) --> 43.9 ms



%%

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '07';...
    '08';...
    '09';...
    '12';...
    '13';...
    '15';...
    '16';...
    ];

for SubjInd = 7:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiFolder = fullfile(SubjectFolder, 'Retinotopy', 'Nifti');
    
    FM_Folder = fullfile(SubjectFolder, 'Retinotopy', 'FieldMap');
    cd(FM_Folder)
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = 43.9;
    
    TEMP = dir('S*_DeltaTE1p02ms_Te7.02.nii');
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase{1} = [fullfile(pwd,TEMP(2).name) ',1'];
    clear TEMP
    
    TEMP = dir('S*_DeltaTE1p02ms_Te6.nii');
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude{1} = [fullfile(pwd,TEMP.name) ',1'];
    clear TEMP
    
    cd(NiftiFolder)
    TEMP=dir; 
    cd(TEMP(3).name); 
    clear TEMP;
    TEMP=dir('S*.nii');
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session.epi{1} = [fullfile(pwd,TEMP(1).name) ',1'];
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = [6 7.02];
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 0;
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = -1;
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = 0;
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = {'home/SHARED/Program/Matlab/spm8/templates/T1.nii'};
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = .5;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = .02;
    
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
    
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
    
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
    
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
    
    
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;
    
    
    cd(FM_Folder)
    
    
    save(strcat('Create_VDM_Subj_', SubjID, '_matlabbatch'));
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end