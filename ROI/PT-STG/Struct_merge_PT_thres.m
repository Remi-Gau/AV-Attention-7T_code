clear; clc

% Folders definitions
RootFolder = fullfile(pwd, '..','..', '..');

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


for SubjInd = 11 % 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    cd(fullfile(RootFolder, 'Subjects_Data'))
        
    Hdr = spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg_data_data.nii']);
    PT_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg_data_data.nii']));
    PT_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_rcr_gm_avg_data_data.nii']));
    
    PT_vol = PT_lcr_vol + PT_rcr_vol;
    
    tabulate(PT_vol(:))
    
    PT_surf = zeros(size(PT_vol));
    PT_surf(PT_vol==1)=1;
    Hdr.fname = ['Subj_' SubjID '_PT_surf_thres.nii'];
    spm_write_vol(Hdr, PT_surf);
   
    cd (RootFolder)
    
end

