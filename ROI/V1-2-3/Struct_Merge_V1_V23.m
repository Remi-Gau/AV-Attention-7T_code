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

    
    cd(fullfile(RootFolder, 'Subjects_Data', 'ROI_V1'))
    
    Hdr = spm_vol(['T1_' SubjID  '_lcr_ProbRet_V1_V2-3_surf_data.nii']);
    V1_V23_lcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_lcr_ProbRet_V1_V2-3_surf_data.nii']));
    V1_V23_rcr_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_rcr_ProbRet_V1_V2-3_surf_data.nii'])); 
    
    V1_V23_vol = V1_V23_rcr_vol + V1_V23_lcr_vol;
    
    unique(V1_V23_vol)
    
    V1_surf = zeros(size(V1_V23_vol));
    V1_surf(V1_V23_vol==22)=1;
    V1_surf(V1_V23_vol==48)=1;
    V1_surf(V1_V23_vol==122)=1;
    V1_surf(V1_V23_vol==144)=1;
    V1_surf(V1_V23_vol==148)=1;
    V1_surf(V1_V23_vol==248)=1;
    Hdr.fname = ['Subj_' SubjID '_V1_surf.nii'];
    spm_write_vol(Hdr, V1_surf);
    
    
    V23_surf = zeros(size(V1_V23_vol));
    V23_surf(V1_V23_vol==26)=1;
    V23_surf(V1_V23_vol==126)=1;
    V23_surf(V1_V23_vol==152)=1;
    Hdr.fname = ['Subj_' SubjID '_V2-3_surf.nii'];
    spm_write_vol(Hdr, V23_surf);


    cd (RootFolder)
    
end

