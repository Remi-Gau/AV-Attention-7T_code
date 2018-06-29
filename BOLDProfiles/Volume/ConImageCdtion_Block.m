%%
clc; clear;

ANTs = 0;

StartFolder=fullfile(pwd, '..', '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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
    
    %     GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'FFX_Block');
    
    %     BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
    %         ['Subject_' SubjID],'UpsampledBetas');
    BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'Transfer');
    

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

    RegNumbers = GetRegNb(SPM);
    
    clear SPM BetaImg
    
    for CondInd = 1:length(Conditions_Names)
        
        tmp=BetaNames(BetaOfInterest,1:length(Conditions_Names{CondInd}));
        
        BetaImgInd{CondInd}=BetaOfInterest(strcmp(Conditions_Names{CondInd}, cellstr(tmp)));  %#ok<*SAGROW>
        
        for i=1:length(BetaImgInd{CondInd})
            if ANTs
                BetaImg{CondInd}(i,:)=fullfile(BackUpFolder, ['rbeta_' sprintf('%3.4d', BetaImgInd{CondInd}(i)) '_def.nii']);
            else
                BetaImg{CondInd}(i,:)=fullfile(BackUpFolder, ['r4beta_' sprintf('%3.4d', BetaImgInd{CondInd}(i)) '.nii']); %#ok<*UNRCH>
            end
        end
        
        H = spm_vol(BetaImg{CondInd});
        V = spm_read_vols(H);
        
        V = nanmean(V,4);
        
        H = H(1);
        H.fname = fullfile(BackUpFolder, ['subj_' SubjID '_mean_' strrep(Conditions_Names{CondInd}, ' ', '_') '.nii']);
        
        spm_write_vol(H,V)
        
        clear tmp i
    end
    clear BetaNames BetaOfInterest
    
    
    
    
end


cd(StartFolder)
