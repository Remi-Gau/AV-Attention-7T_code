%%
clc; clear;

StartFolder=fullfile(pwd, '..', '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

for SubjInd = 2:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');

    BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'UpsampledBetas');
    
    
    %% Opens beta images
    fprintf(' Identifying the relevant beta images\n')
    
    load(fullfile(GLM_Folder, 'SPM.mat'));
    
    BetaNames = char(SPM.xX.name');
    
    BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
        strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
        2);
    
    BetaOfInterest = find(BetaOfInterest);
    
    BetaNames(:,1:6)=[];
    
    clear SPM
    
    for CondInd = 1:length(Conditions_Names)
        
        tmp=BetaNames(BetaOfInterest,1:length(Conditions_Names{CondInd}));
        
        BetaImgInd{CondInd}=BetaOfInterest(strcmp(Conditions_Names{CondInd}, cellstr(tmp)));  %#ok<*SAGROW>
        
        for i=1:length(BetaImgInd{CondInd})
            BetaImg{CondInd}(i,:)=fullfile(BackUpFolder, ['rbeta_' sprintf('%3.4d', BetaImgInd{CondInd}(i)) '.nii']);
        end
        
        clear tmp i
    end
    clear BetaNames BetaOfInterest
    
    
    
    %% AVERAGING
    for CondInd = 1:length(Conditions_Names)/2 % For each Condition
        
        SourceVolumeStruc = spm_vol(BetaImg{CondInd});
        SourceVolume = mean(spm_read_vols(SourceVolumeStruc),4);
        
        SourceVolumeStruc2 = spm_vol(BetaImg{CondInd+3});
        SourceVolume2 = mean(spm_read_vols(SourceVolumeStruc2),4);
        
        FinalVolume = SourceVolume2 - SourceVolume;
        
        hdr = SourceVolumeStruc(1);
        hdr.fname = fullfile(SubjectFolder, ...
            ['Subject_' SubjID '_' strrep(Conditions_Names{CondInd}(1:end-21),' ','_') '-AttMod.nii']); 
        
        spm_write_vol(hdr,FinalVolume)
        
    end
    
    SourceVolumeStruc = spm_vol(cat(1,BetaImg{1},BetaImg{2},BetaImg{3}));
    SourceVolume = mean(spm_read_vols(SourceVolumeStruc),4);
    
    SourceVolumeStruc2 = spm_vol(cat(1,BetaImg{4},BetaImg{5},BetaImg{6}));
    SourceVolume2 = mean(spm_read_vols(SourceVolumeStruc2),4);
    
    FinalVolume = SourceVolume2 - SourceVolume;
    
    hdr = SourceVolumeStruc(1);
    hdr.fname = fullfile(SubjectFolder, ...
        ['Subject_' SubjID '_AllStim-AttMod.nii']);
    
    spm_write_vol(hdr,FinalVolume)    
    
    clear CondInd

    
end


cd(StartFolder)
