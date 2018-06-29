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
    '15';...
    '16'
    ];


for SubjInd = size(SubjectList,1)-1
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);

    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_V2-3'))
    
    Hdr = spm_vol(['T1_' SubjID  '_lcr_inf_V2-3_data.nii']);
    V2_3_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_lcr_inf_V2-3_data.nii']));
    V2_3_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_rcr_inf_V2-3_data.nii'])); 
    
    V2_3_vol = zeros(size(V2_3_lcr_vol));
    
    unique(V2_3_lcr_vol)
    unique(V2_3_rcr_vol)
    
    V2_3_vol(V2_3_lcr_vol==22)=1;
    V2_3_vol(V2_3_lcr_vol==122)=1;
    V2_3_vol(V2_3_rcr_vol==22)=1;
    V2_3_vol(V2_3_rcr_vol==122)=1;
    Hdr.fname = ['Subj_' SubjID '_V2-3_surf.nii'];
    spm_write_vol(Hdr, V2_3_vol);
    
    cd (RootFolder)
    
end

