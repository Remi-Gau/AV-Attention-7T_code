%%
clear
clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

% StartFolder = fullfile(pwd, '..', '..');
StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

DateFormat = 'yyyy_mm_dd_HH_MM';

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


%% Define folders, number of runs, scans per run...
for SubjInd = 1:size(SubjectList,1)
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(SubjectFolder, 'FFX_bu', 'Structural', 'mUNI.nii')};
    
    NbRuns = length( FoldersNames{SubjInd});
    
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {};
    for RunInd=1:NbRuns
        
        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        
        if RunInd==1
            TEMP = dir('meanuS*.nii');
            matlabbatch{1}.spm.spatial.coreg.estimate.source = {fullfile(NiftiSourceFolder, '01', TEMP.name)};
        end
        
        TEMP = dir('uS*.nii');
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1}.spm.spatial.coreg.estimate.other{end+1,1} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        
        cd ..
    end

    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2 1 0.7];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [repmat(0.001, 1, 3), repmat(0.0005, 1, 3), repmat(0.005, 1, 3), repmat(0.0005, 1, 3)];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    cd(fullfile(SubjectFolder, 'FFX_bu', 'Structural'))
    save (strcat('Coreg_RealignedEPI_Subj_', SubjID, '_', datestr(now, DateFormat), 'matlabbatch'));
    spm_jobman('run',matlabbatch)

    
end

cd(StartFolder)

spm_check_registration