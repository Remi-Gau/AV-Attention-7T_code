% compute tSNR image for each fMRI run

% might need double checking as the TPMs might not be registered with the
% raw BOLD though they have same dimension so results should still be OK.

clear; close all; clc;

code_folder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile('/home/remi/github/spmup')))

derivatives_folder = '/home/remi/Dropbox/BIDS/AV_att/derivatives/';
raw_folder = '/home/remi/Dropbox/BIDS/AV_att/rawdata';

output_folder = fullfile(derivatives_folder, 'smpup');

SubLs = dir(fullfile(derivatives_folder, 'spm12', 'sub-*'));
NbSub = numel(SubLs);

% flags for reslicing TPMs
flags.which = 1;
flags.wrap = [1 1 0];
flags.mean = 0;

for iSub = 1:NbSub % for each subject
    
    fprintf(' doing %s\n', SubLs(iSub).name)
    
    spm_mkdir(output_folder, SubLs(iSub).name, {'ses-1' 'ses-2'}, 'func')
    
    output_folder = fullfile(derivatives_folder, 'smpup');
    
    sub_dir = fullfile(derivatives_folder, 'spm12', SubLs(iSub).name);
    
    mean = spm_select('FPList', ...
        fullfile(sub_dir, 'ses-1', 'func'), '^mean.*.nii$');
    
    masks = spm_select('FPList', ...
        fullfile(sub_dir, 'ses-1', 'anat'), '^rc[1-3].*.nii$');
    
    % if the subject TPMs
    if isempty(masks)
        masks = spm_select('FPList', ...
            fullfile(sub_dir, 'ses-1', 'anat'), '^c[1-3].*.nii$');
        
        P = cat(1, cellstr(mean), cellstr(masks));
        
        spm_reslice(P,flags)
    end
    
    % get the bold files
    raw_sub_dir = fullfile(raw_folder, SubLs(iSub).name);
    
    time_series = spm_select('FPListRec', ...
        raw_sub_dir, '^*audiovisualattention.*bold.nii$');
    
    % unzip if necessary
    if size(time_series,1)<4
        time_series = spm_select('FPListRec', ...
            raw_sub_dir, '^*audiovisualattention.*bold.nii.gz$');
        
        for i_run = 1:size(time_series,1)
            gunzip(time_series(i_run,:))
        end
        
        time_series = spm_select('FPListRec', ...
            raw_sub_dir, '^*audiovisualattention.*bold.nii$');
        
    end
    
    for i_func = 1:size(time_series,1)
        
        tSNR = spmup_temporalSNR(time_series(i_func,:), masks, 1);
        
        clear tSNR
        
        [path, file, ext] = spm_fileparts(time_series(i_func,:));
        
        if i_func<3
            ses = '1';
        else
            ses = '2';
        end
        
        movefile(fullfile(pwd, 'tSNR.ps'), ...
            fullfile(output_folder, SubLs(iSub).name, ['ses-' ses], 'func', ['tSNR_' file '.ps']))
        
        tSNR_img = spm_select('FPListRec', ...
            fullfile(raw_sub_dir), '^tSNR.*audiovisualattention.*bold.nii$');
        
        movefile(tSNR_img, ...
            fullfile(output_folder, SubLs(iSub).name, ['ses-' ses], 'func', ['tSNR_' file '.nii']))
    end
    
end


