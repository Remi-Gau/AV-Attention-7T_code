clear; clc

% Folders definitions
RootFolder = pwd;

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


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    ROI_Folder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    
    cd(ROI_Folder)
    
    STG_hdr = spm_vol(fullfile(ROI_Folder, 'rrwHG_STG_AAL.nii'));
    STG_vol = spm_read_vols(STG_hdr);
    
    TE_vol = spm_read_vols(spm_vol(fullfile(ROI_Folder, 'rrwTE_MNI.nii'))); 
    
    STG_vol(TE_vol==1)=0;
    STG_hdr.fname = 'rrwSTG_AAL.nii';
    spm_write_vol(STG_hdr, STG_vol)
    
    [I]=find(TE_vol);
    [X,Y,Z]=ind2sub(STG_hdr.dim,I);
    STG_vol(:,1:min(Y)-1,:)=0;
    STG_hdr.fname = 'rrwSTG_Post_AAL.nii';
    spm_write_vol(STG_hdr, STG_vol)
 
    
    cd (RootFolder)
    
end

