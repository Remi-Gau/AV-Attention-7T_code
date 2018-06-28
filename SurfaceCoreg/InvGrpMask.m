clear; clc

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)

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
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Transfer');
    
    try
        gunzip(fullfile(SubjectFolder, 'mask_bin_*h_data_clone_transform_math_uncrop_def.nii.gz'))
    catch
    end
    
    Masks = spm_select('FPList', fullfile(SubjectFolder),...
            '^mask_bin_.*h_data_clone_transform_math_uncrop_def.nii$');

    hdr = spm_vol(Masks);
    
    vol = spm_read_vols(hdr);
    vol = logical(nansum(ceil(vol),4));
    
    hdr = hdr(1);
    hdr.fname = fullfile(SubjectFolder, ['Subject_' SubjID '_grp_mask_inv.nii']);

    spm_write_vol(hdr,vol)

    clear SubjectFolder SubjID hdr vol Masks
    
end






