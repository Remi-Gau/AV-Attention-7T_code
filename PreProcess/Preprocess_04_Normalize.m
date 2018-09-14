function Preprocess_04_Normalize
% Normalize EPI images
% This is moslty to run a traditional group level analysis to have a look
% at the activations

clear 
clc


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

for SubjInd = 1:size(SubjectList,1)
    
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    
    NbRuns = length(FoldersNames{SubjInd});
    
    
    matlabbatch = {};
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def{1} = fullfile(SubjectFolder, 'FFX', 'Structural', 'y_UNI.nii');
    
    % List the images for each run and the mean image 
    for RunInd=1:NbRuns
        
        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        
        if RunInd==1
            TEMP = dir('mean*.nii');
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = fullfile(NiftiSourceFolder, '01', TEMP.name);
        end
       
        TEMP = dir('URS*.nii');
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{end+1,1} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        
        
    end
    
    char(matlabbatch{1}.spm.spatial.normalise.write.subj.resample)
    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;...
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    cd(fullfile(SubjectFolder, 'FFX', 'Structural'))
    SaveBatch('Normalize_EPI_Subj_',SubjID, matlabbatch)
    
    
    fprintf('\n')
    disp('%%%%%%%%%%%%%%%')
    disp('   NORMALIZE   ')
    disp('%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fullfile(SubjectFolder, 'FFX', 'Structural', 'mUNI.nii')};
    SaveBatch('Normalize_Struc_Subj_',SubjID, matlabbatch)
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end

end

function SaveBatch(Prefix,SubjID, matlabbatch)
DateFormat = 'yyyy_mm_dd_HH_MM';
save (strcat(Prefix, SubjID, '_', datestr(now, DateFormat), '_matlabbatch'), 'matlabbatch');
end

