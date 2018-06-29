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


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);

    
    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_TEs_cyt'))
    
    Hdr = spm_vol(['T1_' SubjID  '_lcr_TE_SubDiv_data.nii']);
    TE_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_lcr_TE_SubDiv_data.nii']));
    TE_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_rcr_TE_SubDiv_data.nii'])); 
    
    TE_vol = TE_rcr_vol + TE_lcr_vol;
    
    unique(TE_vol)
    
    TE12_vol = zeros(size(TE_vol));
    TE12_vol(TE_vol==4)=1;
    Hdr.fname = ['Subj_' SubjID '_TE12.nii'];
    spm_write_vol(Hdr, TE12_vol);
    
    
    TE10_vol = zeros(size(TE_vol));
    TE10_vol(TE_vol==5)=1;
    Hdr.fname = ['Subj_' SubjID '_TE10.nii'];
    spm_write_vol(Hdr, TE10_vol);
    
    
    TE11_vol = zeros(size(TE_vol));
    TE11_vol(TE_vol==6)=1;
    Hdr.fname = ['Subj_' SubjID '_TE11.nii'];
    spm_write_vol(Hdr, TE11_vol);
    
    
%     Hdr.fname = ['Subj_' SubjID '_TE.nii'];
%     spm_write_vol(Hdr, TE_vol>0);

    cd (RootFolder)
    
end

