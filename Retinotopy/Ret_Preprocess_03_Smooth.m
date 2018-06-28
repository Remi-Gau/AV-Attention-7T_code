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
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '15';...
    '16'
    ];

FoldersNames = {...
    1;...
    1:2;...
    1:2;...
    1;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    };

for SubjInd =size(SubjectList,1)-1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Retinotopy', 'Nifti');
    
    NbRuns = length( FoldersNames{SubjInd});
    
    
    matlabbatch = {};
    
    matlabbatch{1}.spm.spatial.smooth.data = {};
    for RunInd=1:NbRuns
        
        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        TEMP = dir('URS*.nii');
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1}.spm.spatial.smooth.data{end+1} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        
    end
    
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'S';
    
    cd(NiftiSourceFolder)
    save (strcat('SmoothNonNormed_Subj_', SubjID, '_matlabbatch'));
    
    
    fprintf('\n')
    disp('%%%%%%%%%%%%')
    disp('   SMOOTH   ')
    disp('%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end

% matlabpool('CLOSE')
