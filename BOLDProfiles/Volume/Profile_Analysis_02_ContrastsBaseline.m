clear; clc;

% StartDirectory = fullfile(pwd, '..','..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';
cd (StartDirectory)

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

CondNames = {...
    'A Stim - Auditory Attention';....
    'V Stim - Auditory Attention';....
    'AV Stim - Auditory Attention';....
    'A Stim - Visual Attention';....
    'V Stim - Visual Attention';....
    'AV Stim - Visual Attention';....
    'Baseline'}


for SubjInd = 2:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    load(fullfile(SubjectFolder, 'FFX_baseline', 'SPM.mat'))
    BetaNames = char(SPM.xX.name');
    
    

    AnalysisFolder = fullfile(SubjectFolder, 'FFX_baseline', 'Beta');
    
    cd(AnalysisFolder)
    
    for iCond = 1%:numel(CondNames)
        tmp = strfind(cellstr(BetaNames),[' ' CondNames{iCond} '*bf(1)']);
        tmp = ~cellfun('isempty',tmp);
        BetaOfInterest = find(tmp);
        for iIMG = 1:numel(BetaOfInterest)
            IMGs(iIMG,:) = spm_select('FPList', pwd,...
                ['^rbeta_' sprintf('%04.0f',BetaOfInterest(iIMG))  '.nii$']);
        end
        
        IMGs
        
        hdr = spm_vol(IMGs);
        vol = spm_read_vols(hdr);
         
        hdr=hdr(1);
        hdr.fname=fullfile(pwd,[strrep(CondNames{iCond},' ','') '.nii']);
        
        spm_write_vol(hdr,mean(vol,4));
        
    end

end