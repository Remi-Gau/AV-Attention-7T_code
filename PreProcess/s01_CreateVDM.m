% Creates a voxel displacement map starting with the fieldmap

clear 
clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

% Folders definitions
% RootFolder = fullfile(pwd, '..', '..');
RootFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

%% Phase enconded lines (PELines) and ReadOutTime
% PELines = ((BaseResolution * PartialFourier)/iPat) + ((iPat-1)/iPAT) * ReferenceLines) =
% ReadoutDuration = PELines * InterEchoSpacing

% GRAPPA=iPAT4 ; Partial Fourrier=6/8 ; 48 sli ; TE=25ms ; Res=0.75 mm
% Bandwidth Per Pixel Phase Encode = 15.873


%% According to R. Trampel

% For distortion correction: ignore Partial Fourrier and references lines
% BaseResolution/iPAT = PELines

% Effective echo spacing: 2 ways to calculate, should be the same
% 1/(Bandwidth Per Pixel Phase Encode * Reconstructed phase lines) -->  0.246 ms
% echo spacing (syngo) / iPAT

% SPM Total readout time = 1/"Bandwidth Per Pixel Phase Encode", stored in
% DICOM tag (0019, 1028) --> 63 ms



%%

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

% List the number of folders for each subject.
% Som subject's had fewer sessions
FoldersNames = {...
    1:4;...
    1:4;...
    1:4;...
    1:2;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:2;...
    1:4;...
    1:4;...
    };

% List the field map folders
FM_FolderList = {'01';'02'};


for SubjInd = 1:size(SubjectList,1)
    
    matlabbatch = {};
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    TargetFile = dir(fullfile(NiftiFolder,'01','S*_iPAT4_6_8_48sli_TE25_0p75_Te25.nii'));
    TargetFile = fullfile(NiftiFolder,'01',TargetFile.name);
    
    FM_Folder = fullfile(SubjectFolder, 'FieldMap');
    
    NbFM = numel(FoldersNames{SubjInd})/2;
    
    % for each field map
    for FM_Ind = 1:NbFM
        
        cd(fullfile(FM_Folder,FM_FolderList{FM_Ind}))
        
        %% Coregister to the first image of the first run  
        % This is required as the VDM will be applied after realignement
        % that gets the very first image as reference
        TEMP = dir('S*_DeltaTE1p02ms_Te6.nii');        
        matlabbatch{end+1}.spm.spatial.coreg.estimate.source{1} = ...
            fullfile(FM_Folder,FM_FolderList{FM_Ind},[TEMP.name ',1']);
        
        matlabbatch{end}.spm.spatial.coreg.estimate.ref{1} = [TargetFile ',1'];
        
        TEMP = dir('S*_DeltaTE1p02ms_Te7.02.nii');
        matlabbatch{end}.spm.spatial.coreg.estimate.other{1} = ...
            fullfile(FM_Folder,FM_FolderList{FM_Ind},[TEMP.name ',1']);
        
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [8 4 2 1];
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [repmat(0.001, 1, 3), repmat(0.0005, 1, 3), repmat(0.005, 1, 3), repmat(0.0005, 1, 3)];
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

        
        %% Batch to create VDM
        matlabbatch{end+1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = 63; %#ok<*SAGROW>
        
        TEMP = dir('S*_DeltaTE1p02ms_Te7.02.nii');
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.phase{1} = ...
            [fullfile(FM_Folder,FM_FolderList{FM_Ind},TEMP.name) ',1'];
        clear TEMP
        
        TEMP = dir('S*_DeltaTE1p02ms_Te6.nii');
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.magnitude{1} = ...
            [fullfile(FM_Folder,FM_FolderList{FM_Ind},TEMP.name) ',1'];
        clear TEMP

        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.session.epi{1} = [TargetFile ',1'];
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = [6 7.02];
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 0;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = -1;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = 0;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = 'Mark3D';
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = ...
            {'home/SHARED/Program/Matlab/spm8/templates/T1.nii'};
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = .5;
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = .02;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 0;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
        
        matlabbatch{end}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;

        cd(FM_Folder)
        
    end
    
    save(strcat('Create_VDM_Subj_', SubjID, '_matlabbatch'));
    
    spm_jobman('run',matlabbatch)

    cd (RootFolder)
    
end