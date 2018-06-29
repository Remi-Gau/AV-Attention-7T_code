clear; clc

% Folders definitions
RootFolder = fullfile(pwd, '..','..','..');

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


for SubjInd = [4 11] %size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);

    
    cd(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1'))
    
    Hdr = spm_vol(['Subj_' SubjID  '_A1_lcr_norm_RG_UN_data.nii']);
    A1_lcr_vol = spm_read_vols(spm_vol(['Subj_' SubjID  '_A1_lcr_norm_RG_UN_data.nii']));
    A1_rcr_vol = spm_read_vols(spm_vol(['Subj_' SubjID  '_A1_rcr_norm_RG_UN_data.nii'])); 
    
    A1_vol = zeros(size(A1_lcr_vol));
    
    unique(A1_lcr_vol)
    unique(A1_rcr_vol)
    
    A1_vol(A1_lcr_vol>0)=1;
    A1_vol(A1_rcr_vol>0)=1;
    Hdr.fname = ['T1_' SubjID '_A1_surf.nii'];
    spm_write_vol(Hdr, A1_vol);
    
    

    cd (RootFolder)
    
end

