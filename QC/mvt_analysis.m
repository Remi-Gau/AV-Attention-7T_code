% Analyizes realignement parameters and computes / plots several QC results
clear; clc; close

code_folder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile('/home/remi/github/spmup')))

data_folder = '/media/remi/BackUp2/AV_Integration_7T_2/derivatives/';

output_folder = fullfile(data_folder, 'smpup');
mkdir(output_folder)
output_mat = fullfile(output_folder, 'mvt_recap.mat');
output_tsv = fullfile(output_folder, 'mvt_recap.csv');

SubLs = dir(fullfile(data_folder, 'spm12', 'sub-*'));
NbSub = numel(SubLs);

FD_threshold = .4; % mm

idx = 1;

for iSub = 1:NbSub % for each subject
    
    fprintf('Processing %s\n', SubLs(iSub).name)
    
    % Subject directory
    SubDir = fullfile(data_folder, 'spm12', SubLs(iSub).name);
    Anat_spm_dir = fullfile(SubDir, 'ses-1', 'anat');
    
    % Compute radius of brain center to surface
    anat_file = spm_select('FPList', Anat_spm_dir, ['^' SubLs(iSub).name  '_ses-1_T1w.nii$']);
    radius = spmup_comp_dist2surf(anat_file);
    
    % Plots results
    RP_files = spm_select('FPListRec', SubDir , '^rp.*.txt$');
    for iRP_file = 1:size(RP_files,1)
        [FD,RMS,motion] = spmup_FD(RP_files(iRP_file,:),radius);
        
        MovementRecap(idx,1:2) = [iSub iRP_file];  %#ok<*SAGROW>
        % mean framewise displacement for this run
        MovementRecap(idx,3) = mean(FD); 
        % number of time points with FD above threshold
        MovementRecap(idx,4) = sum(FD>FD_threshold); 
        % proportion of time points with FD above threshold
        MovementRecap(idx,5) = sum(FD>FD_threshold)/size(FD,1); 
        % mean root mean square of mvt 
        MovementRecap(idx,6) = mean(RMS);
        % number of outlier time points for RMS
        MovementRecap(idx,7) = sum(spmup_comp_robust_outliers(RMS)); 
        % proportion of outlier time points for RMS
        MovementRecap(idx,8) = sum(spmup_comp_robust_outliers(RMS))/size(RMS,1); 
        
        idx =  idx +1;
        
        [path, file, ext] = spm_fileparts(RP_files(iRP_file,:));
        
        movefile(fullfile(path, 'displacement.ps'), ...
            fullfile(output_folder, [file '.ps']))
    end
    
end




%%
Legends = {...
     'subject - run', ...
     'mean FD', 'FD outliers #time points' 'FD outliers %time points',...
     'mean RMS', 'RMS outliers #time points' 'RMS outliers %time points'};

save(output_mat, 'Legends', 'MovementRecap', 'SubLs')

fid = fopen (output_tsv, 'w');

for i=1:length(Legends)
    fprintf (fid, '%s,', Legends{i});
end

for iRow=1:size(MovementRecap,1)
    fprintf (fid, '\n');
    fprintf (fid, [SubLs(MovementRecap(iRow,1)).name '_run-' num2str(MovementRecap(iRow,2)) ',']);
    for iCol=3:size(MovementRecap,2)
        fprintf (fid, '%f,', MovementRecap(iRow,iCol));  
    end
end


fprintf (fid, '\n\n\n\n');

Legends = {...
     'subject', ...
     'mean FD', 'FD outliers #time points' 'FD outliers %time points',...
     'mean RMS', 'RMS outliers #time points' 'RMS outliers %time points'};

 for i=1:length(Legends)
    fprintf (fid, '%s,', Legends{i});
end

 for iSub = 1:NbSub
    fprintf (fid, '\n');
    fprintf (fid, [SubLs(iSub).name ',']);
    row_to_select = MovementRecap(:,1)==iSub;
    for iCol=3:size(MovementRecap,2)
        MovementRecapSubj(iSub,iCol-2) =  mean(MovementRecap(row_to_select,iCol));
        fprintf (fid, '%f,', mean(MovementRecap(row_to_select,iCol)));  
    end
 end

fprintf (fid, '\n\n\nMEAN,'); 
for iCol=1:size(MovementRecapSubj,2)
    fprintf (fid, '%f,', mean(MovementRecapSubj(:,iCol)));
end
fprintf (fid, '\nSTD,'); 
for iCol=1:size(MovementRecapSubj,2)
    fprintf (fid, '%f,', std(MovementRecapSubj(:,iCol)));
end

fclose (fid);