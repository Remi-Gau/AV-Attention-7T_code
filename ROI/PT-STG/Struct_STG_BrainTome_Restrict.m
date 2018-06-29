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
    
    cd(fullfile(SubjectFolder, 'Transfer', 'ROI'))
    
    TE_vol = spm_read_vols(spm_vol(['T1_' SubjID  '_A1_surf.nii']));

    ImgLs = dir('rrw*Tome.nii');
    
    for iImg=1:numel(ImgLs)
        
        STG_hdr = spm_vol(ImgLs(iImg).name);
        STG_vol = spm_read_vols(STG_hdr);
        
        STG_vol(TE_vol==1)=0;
        
        [I]=find(TE_vol);
        [X,Y,~]=ind2sub(STG_hdr.dim,I);
        Y(X<round(size(STG_vol,1)/2,1))=[];
        STG_vol(1:round(size(STG_vol,1)/2,1), 1:min(Y)-1, :) = 0;
        
        [I]=find(TE_vol);
        [X,Y,~]=ind2sub(STG_hdr.dim,I);
        Y(X>round(size(STG_vol,1)/2,1))=[];
        STG_vol(round(size(STG_vol,1)/2,1):end, 1:min(Y)-1, :) = 0;
        
        STG_hdr.fname = ['rrw_p' STG_hdr.fname(4:end)];
            spm_write_vol(STG_hdr, STG_vol);
        
    end
    
    cd (RootFolder)
    
end

