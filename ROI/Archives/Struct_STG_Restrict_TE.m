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

    
    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_TE'))
    TE_lcr_vol = spm_read_vols(spm_vol(['Subj_' SubjID  '_lcr_TE_RG_data.nii']));
    TE_rcr_vol = spm_read_vols(spm_vol(['Subj_' SubjID  '_rcr_TE_RG_data.nii'])); 
    
    
    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_STG'))
    STG_hdr = spm_vol(['Subj_' SubjID '_HG_STG_AAL.nii']);
    STG_vol = spm_read_vols(STG_hdr);
    
        
    STG_vol(TE_lcr_vol==1)=0;
    STG_vol(TE_rcr_vol==1)=0;
    STG_hdr.fname = ['Subj_' SubjID '_STG_AAL.nii'];
    spm_write_vol(STG_hdr, STG_vol);
    
    
    [I]=find(TE_lcr_vol);
    [X,Y,Z]=ind2sub(STG_hdr.dim,I);
    STG_vol(1:round(size(STG_vol,1),1), 1:min(Y)-1, :) = 0;
    
    [I]=find(TE_rcr_vol);
    [X,Y,Z]=ind2sub(STG_hdr.dim,I);
    STG_vol(round(size(STG_vol,1),1):end, 1:min(Y)-1, :) = 0;
    
    STG_hdr.fname = ['Subj_' SubjID '_pSTG.nii'];
    spm_write_vol(STG_hdr, STG_vol);
 
    
    cd (RootFolder)
    
end

