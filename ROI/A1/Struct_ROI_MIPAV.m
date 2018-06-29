%%
clear; clc;

StartFolder = pwd;

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

%%
for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    fprintf(['\n\nProcessing subject : ' SubjID '\n'])

%     mkdir(fullfile(SubjectFolder, 'ROI_MIPAV'))
%     
%     cd(fullfile(StartFolder, 'Subjects_Data'))
%     
%     movefile(['T1_' SubjID '*data.nii'], fullfile(SubjectFolder, 'ROI_MIPAV'))
    
    cd (fullfile(SubjectFolder, 'ROI_MIPAV'))
      
    fprintf('\n Creating ROI')
    
    % Left
    fname = dir(['T1_' SubjID '_lcr*data.nii']);
    
    hdr = spm_vol(fname.name); clear fname
    
    vol_L = spm_read_vols(hdr);
    vol_L(vol_L~=35)=0;
    vol_L(vol_L==35)=1;
    vol_L=logical(vol_L);
    
    hdr = struct(...
        'fname',   'A1_L.nii',...
        'dim',     size(vol_L),...
        'dt',      [spm_type('uint8') spm_platform('bigend')],...
        'mat',     hdr.mat,...
        'pinfo',   [1 0 0]',...
        'descrip', 'A1_L');
    
    hdr = spm_data_hdr_write(hdr);
    
    spm_write_vol(hdr,vol_L); clear hdr vol;

    % Right
    fname = dir(['T1_' SubjID '_rcr*data.nii']);
    
    hdr = spm_vol(fname.name); clear fname
    
    vol_R = spm_read_vols(hdr);
    vol_R(vol_R~=135)=0;
    vol_R(vol_R==135)=1;
    vol_R=logical(vol_R);
    
    hdr = struct(...
        'fname',   'A1_R.nii',...
        'dim',     size(vol_R),...
        'dt',      [spm_type('uint8') spm_platform('bigend')],...
        'mat',     hdr.mat,...
        'pinfo',   [1 0 0]',...
        'descrip', 'A1_R');
    
    hdr = spm_data_hdr_write(hdr);
    
    spm_write_vol(hdr,vol_R);
    

    % Both
    hdr = struct(...
        'fname',   'A1.nii',...
        'dim',     size(vol_R),...
        'dt',      [spm_type('uint8') spm_platform('bigend')],...
        'mat',     hdr.mat,...
        'pinfo',   [1 0 0]',...
        'descrip', 'A1');
    
    hdr = spm_data_hdr_write(hdr);
    vol = sum(cat(4,vol_R,vol_L),4);
    spm_write_vol(hdr,vol);
    
    
    % Both only with cortical voxels
    Cx = logical(spm_read_vols(spm_vol(fullfile(fullfile(SubjectFolder,'Structural', ...
        'CBS', 'Layering', 'T1_04_Layers.nii')))));
    hdr.fname = 'A1_Cx.nii';
    vol(Cx==0)=0;
    spm_write_vol(hdr,vol); 
    clear hdr Cx vol
    
end


fprintf('\n\n')

cd(StartFolder)
