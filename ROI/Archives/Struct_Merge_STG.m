clear; clc

% Folders definitions
RootFolder = fullfile(pwd, '..','..');

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
%     '15';...
    '16'
    ];


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);

    
    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_STG'))
    
    Hdr = spm_vol(['T1_' SubjID  '_lcr_pSTG_data.nii']);
    STG_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_lcr_pSTG_data.nii']));
    STG_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_rcr_pSTG_data.nii'])); 
    
    STG_vol = zeros(size(STG_lcr_vol));
    
    STG_vol(STG_lcr_vol==31)=1;
    STG_vol(STG_lcr_vol==131)=1;
    STG_vol(STG_rcr_vol==31)=1;
    STG_vol(STG_rcr_vol==131)=1;
    Hdr.fname = ['Subj_' SubjID '_STG.nii'];
    spm_write_vol(Hdr, STG_vol);
    
    

    cd (RootFolder)
    
end

