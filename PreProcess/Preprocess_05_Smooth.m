function Preprocess_05_Smooth
% Smooths the EPI images
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

for SubjInd = 1:size(SubjectList,1)
    
    
    SubjID = SubjectList(SubjInd,:)
    
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    
    NbRuns = length(FoldersNames{SubjInd});
    
    
    matlabbatch = {};
    
    matlabbatch{1}.spm.spatial.smooth.data = {};
    % List the images for each run
    for RunInd=1:NbRuns
        
        RunFolder =  ...
            fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd)));
        
        TEMP = dir(fullfile(RunFolder,'wURS*.nii'));
        TEMP = spm_vol(fullfile(RunFolder, TEMP.name));
        
        for i=1:length(TEMP)
            matlabbatch{1}.spm.spatial.smooth.data{end+1,1} = ...
                [TEMP(i).fname ,',', num2str(i)];
        end
        
    end
    
    matlabbatch{1}.spm.spatial.smooth.fwhm = [2.5 2.5 2.5];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    SaveBatch(fullfile(SubjectFolder,'Smooth_EPI_Subj_'), SubjID, matlabbatch)
    
    
    fprintf('\n')
    disp('%%%%%%%%%%%%')
    disp('   SMOOTH   ')
    disp('%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end


end


function SaveBatch(Prefix, SubjID, matlabbatch)
DateFormat = 'yyyy_mm_dd_HH_MM';
save (strcat(Prefix, SubjID, '_', datestr(now, DateFormat), '_matlabbatch.mat'), 'matlabbatch');
end
