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


for SubjInd = 3 %1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    cd(fullfile(RootFolder, 'Subjects_Data'))
        
    Hdr = spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg_data_data.nii']);
    V1_V23_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg_data_data.nii']));
    V1_V23_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_thresh_clone_transform_strip_clone_transform_bound_mems_rcr_gm_avg_data_data.nii']));
    
    V1_V23_vol = V1_V23_rcr_vol + V1_V23_lcr_vol;
    
    tabulate(V1_V23_vol(:))
    
    V1_surf = zeros(size(V1_V23_vol));
    V1_surf(V1_V23_vol==1)=1;
    Hdr.fname = ['Subj_' SubjID '_V1_surf_thres.nii'];
    spm_write_vol(Hdr, V1_surf);
    
    
    V23_surf = zeros(size(V1_V23_vol));
    V23_surf(V1_V23_vol==2)=1;
    Hdr.fname = ['Subj_' SubjID '_V2-3_surf_thres.nii'];
    spm_write_vol(Hdr, V23_surf);


    cd (RootFolder)
    
end

