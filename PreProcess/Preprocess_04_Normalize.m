function Preprocess_04_Normalize
% Normalize EPI images
% This is moslty to run a traditional group level analysis to have a look
% at the activations

clear 
clc

%  Folders definitions
RootFolder = '/media/remi/BackUp2/AV_Integration_7T_2';

SubjectList = [...
    '02';...
    '03';...
    '04';...
%     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
%     '14';...
    '15';...
    '16'
    ];

FoldersNames = {...
    1:4;...
    1:4;...
    1:4;...
%     1:2;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
%     1:2;...
    1:4;...
    1:4;...
    };

for SubjInd = 2:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    NbRuns = length(FoldersNames{SubjInd});
    
    matlabbatch = {};
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def{1} = ...
        fullfile(SubjectFolder, 'Structural', 'SPM_Segmentation', 'y_UNI.nii');
    
    % List the images for each run and the mean image 
    for RunInd=1:NbRuns
        
        RunFolder = fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd)));

        if RunInd==1
            TEMP = dir(fullfile(RunFolder,'meanUR*.nii'));
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = ...
                fullfile(RunFolder, TEMP.name);
        end
       
        TEMP = dir(fullfile(RunFolder,'URS*.nii'));
        TEMP = spm_vol(fullfile(RunFolder, TEMP.name));
        
        for i=1:length(TEMP)
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{end+1,1} = ...
                [TEMP(i).fname ,',', num2str(i)];
        end
        
        
    end
    
    char(matlabbatch{1}.spm.spatial.normalise.write.subj.resample)
    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [...
        -78 -112 -70;...
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [.75 .75 .75];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    SaveBatch(fullfile(SubjectFolder,'Normalize_EPI_Subj_'), SubjID, matlabbatch)
    
    
    fprintf('\n')
    disp('%%%%%%%%%%%%%%%')
    disp('   NORMALIZE   ')
    disp('%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
%     matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ...
%         {fullfile(SubjectFolder, 'Structural', 'SPM_Segmentation', 'mUNI.nii')};
%     
%     SaveBatch(fullfile(SubjectFolder, 'Normalize_Struc_Subj_'), SubjID, matlabbatch)
%     
%     spm_jobman('run',matlabbatch)

end

end

function SaveBatch(Prefix, SubjID, matlabbatch)
DateFormat = 'yyyy_mm_dd_HH_MM';
save (strcat(Prefix, SubjID, '_', datestr(now, DateFormat), '_matlabbatch.mat'), 'matlabbatch');
end

