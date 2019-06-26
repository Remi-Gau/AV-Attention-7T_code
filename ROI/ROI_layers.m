% Creates a ROI of the lyr 
clear; clc

% Folders definitions
data_folder = '/media/remi/BackUp2/AV_Integration_7T_2/derivatives';

subs = {...
    '02'
    '03'
    '04'
%     '06'
    '07'
%     '08'
%     '09'
    '11'
    '12'
%     '13'
%     '14'
%     '15'
%     '16'
};


ROIs= {...
    'A1-surf'
    'PT-surf-thres'
    'V1-surf-thres'
    'V2-3-surf-thres'
    };


for i_sub = 1:numel(subs)
    
    roi_folder = fullfile(data_folder, 'cbs_tools', ['sub-' subs{i_sub}], ...
        'roi');
    
    layer_folder = fullfile(data_folder, 'cbs_tools', ['sub-' subs{i_sub}], ...
        'ses-1', 'anat');
    
    layer_file = spm_select('FPList', layer_folder, ....
        '^.*ses-1_T1map_bound_layering-6_mems_cr_gm_labels.nii.gz$');
    
    gunzip(layer_file)
    
    layer_file = spm_select('FPList', layer_folder, ....
        '^.*ses-1_T1map_bound_layering-6_mems_cr_gm_labels.nii$');
    
    Hdr = spm_vol(layer_file);
    Vol = spm_read_vols(Hdr);
    
    for i_roi=1:numel(ROIs)
        
        ROI_img = spm_select('FPList', roi_folder, ['^sub-' subs{i_sub} '_roi-' ROIs{i_roi} '.nii$']);
        
        ROI = logical(spm_read_vols(spm_vol(ROI_img)));
        
        Layers = zeros(size(ROI));
        Layers(ROI) = Vol(ROI);
        
        Hdr.fname = fullfile(roi_folder, ['sub-' subs{i_sub} '_roi-' ROIs{i_roi} '_layers-6.nii']);
        
        spm_write_vol(Hdr, Layers)
    end

end

