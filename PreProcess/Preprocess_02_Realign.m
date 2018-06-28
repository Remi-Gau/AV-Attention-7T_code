function Preprocess_02_Realign
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

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = [];
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');

    NbRuns = length( FoldersNames{SubjInd});
    
    
    %% ---------------------------  %
    %      UNWARP & REALIGN         %
    %  ---------------------------  %
    
    tic
    
    matlabbatch = {};
    
    for RunInd=1:NbRuns
        
        RunInd
        
        TEMP = [];
        
        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        TEMP = dir('S1*iPAT4_6_8_48sli_TE25_0p75_Te25.nii');
%         copyfile(TEMP.name, ['bu_' TEMP.name])
        
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1,1}.spm.spatial.realign.estimate.data{1,RunInd}{i,1} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        
    end
    
    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 2; % Defaut = 4
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 3; % Default = 5
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % Register to mean: 0 ---> FALSE
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 1 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
    

    cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(1))))
    SaveBatch(SubjID, matlabbatch)
    
    fprintf('\n\n')
    disp('%%%%%%%%%%%%%%%%')
    disp('   REALIGNING   ')
    disp('%%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
    toc
    
end

end

function SaveBatch(SubjID, matlabbatch)
save (strcat('Realign_Subj_', SubjID, '_matlabbatch'), 'matlabbatch');
end