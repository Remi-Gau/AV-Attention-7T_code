% Gets average smoothness of raw data
clear; clc; close

code_folder = '/home/remi/github/AV-Attention-7T_code';

data_folder = '/home/remi/Dropbox/BIDS/AV_att/derivatives/';

output_folder = fullfile(data_folder, 'afni');

output_tsv = fullfile(output_folder, 'FWHM_raw.csv');

SubLs = dir(fullfile(data_folder, 'afni', 'sub-*'));
NbSub = numel(SubLs);

FWHM_recap = [];

% read data from each run and average across volumes
for iSub = 1:NbSub
    fprintf('Processing %s\n', SubLs(iSub).name)
    
    % Subject directory
    SubDir = fullfile(data_folder, 'afni', SubLs(iSub).name);
    
    % get result file for each run
    file_ls = spm_select('FPList', SubDir , '^FWHM.*.txt$');
    
    FWHM_subject = [];
    
    for i_file = 1:size(file_ls, 1)
        FWHM = load(file_ls(i_file,:)); % load data
        FWHM_subject(end+1, :) = mean(FWHM);
        clear FWHM
    end
    
    FWHM_recap(end+1, :) = mean(FWHM_subject);
end


%  write output file
fid = fopen (output_tsv, 'w');

for iSub = 1:NbSub
    fprintf (fid, '%f,%f,%f\n', FWHM_recap(iSub,:));
end

fclose (fid);