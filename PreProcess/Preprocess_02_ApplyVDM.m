function Preprocess_02_ApplyVDM
% Only applies the voxel displacement map

clear; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


%  Folders definitions
% RootFolder = fullfile(pwd, '..', '..');
RootFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

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


FieldMaps = {...
    [1 1 2 2];... %02
    [1 1 2 2];... %03
    [1 1 2 2];... %04
    [1 1];... %06
    [1 1 2 2];... %07
    [1 1 2 2];... %08
    [1 1 2 2];... %09
    [1 1 2 2];... %11
    [1 1 2 2];... %12
    [1 1 2 2];... %13
    [1 1];... %14
    [1 1 2 2];... %15
    [1 1 2 2];... %16
    };


for SubjInd = 2:size(SubjectList,1)
    
    SubjID = [];
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    FieldmapsFolder = fullfile(SubjectFolder, 'FieldMap');
    
    NbRuns = length( FoldersNames{SubjInd});
    
    
    %  ---------------  %
    %      APPLY VDM    %
    %  ---------------  %
    
    tic
    
    matlabbatch = {};
    
    for RunInd=1:NbRuns
        
        RunInd

        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        TEMP = dir('S1*iPAT4_6_8_48sli_TE25_0p75_Te25.nii');
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1,1}.spm.tools.fieldmap.applyvdm.data(1,RunInd).scans{i,1} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        TEMP = []; %#ok<*NASGU>
        
        cd(fullfile(FieldmapsFolder, sprintf('%2.2d', FieldMaps{SubjInd}(RunInd))))
        TEMP=dir('vdm*.nii');
        matlabbatch{1,1}.spm.tools.fieldmap.applyvdm.data(1,RunInd).vdmfile = {fullfile(pwd, [TEMP.name ',1'])};
        TEMP = [];
        
    end
    
    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 1 0];
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
    matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'u2';
    
    cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(1))))
    SaveBatch(SubjID, matlabbatch)
    
    fprintf('\n\n')
    disp('%%%%%%%%%%%%%%%')
    disp('   Apply VDM   ')
    disp('%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
    toc
    
end

end

function SaveBatch(SubjID, matlabbatch)
save (strcat('ApplyVDM_Subj_', SubjID, '_matlabbatch'), 'matlabbatch');
end