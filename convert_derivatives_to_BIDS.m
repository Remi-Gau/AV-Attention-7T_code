% small script to convert relevant preprocessed data to almost BIDS derivatives

clear
clc
close all

%     spm_mkdir(target_folder, 'cbs_tools', ['sub-' subs{i_sub}], 'ses-1', 'anat')

source_folder = '/media/remi/BackUp2/AV_Integration_7T_2/Subjects_Data';
target_folder = '/media/remi/BackUp2/AV_Integration_7T_2/derivatives';

subs = {...
    '02'
    '03'
    '04'
    '06'
    '07'
    '08'
    '09'
    '11'
    '12'
    '13'
    '14'
    '15'
    '16'};

ROIs= {...
    'A1_surf'
    'PT_surf_thres'
    'V1_surf_thres'
    'V2-3_surf_thres'
    };

for i_sub = 1:numel(subs)
    
    %% anat SPM
    sub_folder = fullfile(source_folder, ['Subject_' subs{i_sub}], ...
        'Structural', 'SPM_Segmentation');
    
    dest_folder = fullfile(target_folder, 'spm12', ['sub-' subs{i_sub}], ...
        'ses-1', 'anat');
    
    img_ls = spm_select('FPList',sub_folder,'^*.nii$');
    mat_ls = spm_select('FPList',sub_folder,'^*UNI.*.mat$');
    
    pattern = 'UNI';
    replace = ['sub-' subs{i_sub} '_ses-1_T1w'];
    
    % transfer nifti
    for i_img = 1:size(img_ls,1)
        
        SOURCE = img_ls(i_img,:);
        
        [path, file, ext] = spm_fileparts(SOURCE);
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    % transfer .mat
    for i_mat = 1:size(mat_ls,1)
        
        SOURCE = mat_ls(i_mat,:);
        
        [path, file, ext] = spm_fileparts(SOURCE);
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    
    %% anat CBS
    sub_folder = fullfile(source_folder, ['Subject_' subs{i_sub}], ...
        'Structural', 'CBS');
    
    dest_folder = fullfile(target_folder, 'cbs_tools', ['sub-' subs{i_sub}], ...
        'ses-1', 'anat');
    
    img_ls = spm_select('FPList',sub_folder,'^*thres.*.nii$');
    
    vtk_ls = spm_select('FPList',sub_folder,'^*thres.*.vtk$');
    
    layering_ls = spm_select('FPList',...
        fullfile(sub_folder, 'Layering', '6', 'VolumetricLayering'), ...
        '^*thres.*.nii.gz$');
    
    seg_ls = spm_select('FPList',...
        fullfile(sub_folder, 'Segmentation', 'exp*', 'exp*-CADAA', 'MGDMMultiBrainSegmentation'), ...
        '^*thres.*.nii.gz$');
    
    pattern_1 = '_thresh_clone_transform_strip_clone_transform';
    pattern_2 = ['T1_' subs{i_sub}];
    replace = ['sub-' subs{i_sub} '_ses-1_T1map'];
    
    % transfer nifti
    for i_img = 1:size(img_ls,1)
        
        SOURCE = deblank(img_ls(i_img,:));
        
        [path, file, ext] = spm_fileparts(SOURCE);
        file = strrep(file, pattern_1, '');
        file = strrep(file, pattern_2, replace);
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    % transfer vtk
    for i_vtk = 1:size(vtk_ls,1)
        
        SOURCE = deblank(vtk_ls(i_vtk,:));
        
        [path, file, ext] = spm_fileparts(SOURCE);
        file = strrep(file, pattern_1, '');
        file = strrep(file, pattern_2, replace);
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    % transfer layering
    pattern_1 = '_thresh_clone_transform_strip_clone_transform_bound';
    
    for i_layer = 1:size(layering_ls,1)
        
        SOURCE = deblank(layering_ls(i_layer,:));
        
        [path, file, ext] = spm_fileparts(SOURCE);
        file = strrep(file, pattern_1, '_bound_layering-6');
        file = strrep(file, pattern_2, replace);
        file = strrep(file, '__', '_');
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    % transfer segmentation
    pattern_1 = '_thresh_clone_transform_strip_clone_transform';
    
    for i_seg = 1:size(seg_ls,1)
        
        SOURCE = deblank(seg_ls(i_seg,:));
        
        [path, file, ext] = spm_fileparts(SOURCE);
        file = strrep(file, pattern_1, '');
        file = strrep(file, pattern_2, replace);
        file = strrep(file, '__', '_');
        
        DESTINATION = fullfile( dest_folder, ...
            [strrep(file, pattern, replace) ext] );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    end
    
    % transfer ROIs
    sub_folder = fullfile(source_folder, ['Subject_' subs{i_sub}], ...
        'ROI');
    
    dest_folder = fullfile(target_folder, 'cbs_tools', ['sub-' subs{i_sub}], ...
        'roi');
    
    mkdir(dest_folder)
    
    pattern_1 = '_';
    pattern_2 = '-';
    replace = ['sub-' subs{i_sub} '_roi-'];
    
    for i_roi = 1:numel(ROIs)
        
        SOURCE = spm_select('FPList', sub_folder, ['^' ROIs{i_roi} '.nii$']);
        
        [path, file, ext] = spm_fileparts(SOURCE);
        file = strrep(file, pattern_1, pattern_2);
        
        DESTINATION = fullfile( dest_folder, ...
            [replace file ext] );
        
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE, DESTINATION)
    end
    
    

    
    
    %% fieldmaps
    sub_folder = fullfile(source_folder, ['Subject_' subs{i_sub}], ...
        'FieldMap');
    
    dest_folder = fullfile(target_folder, 'spm12', ['sub-' subs{i_sub}]);
    
    vdm_ls = spm_select('FPListRec',...
        fullfile(sub_folder), ...
        '^vdm.*.nii$');
    
    for i_vdm = 1:size(vdm_ls,1)
        
        SOURCE = deblank(vdm_ls(i_vdm,:));
        
        file = ['vdm_scsub-' subs{i_sub} '_ses-' num2str(i_vdm) '_acq-48_run-01_phasediff.nii'];
        
        DESTINATION = fullfile( dest_folder, ['ses-' num2str(i_vdm)], ...
            file );
        
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION)
    end
    
    
    %% func
    sub_folder = fullfile(source_folder, ['Subject_' subs{i_sub}], ...
        'Nifti', 'NoMoCo');
    
    dest_folder = fullfile(target_folder, 'spm12', ['sub-' subs{i_sub}]);
    
    func_ls = spm_select('FPListRec',...
        fullfile(sub_folder), ...
        '^URS.*.nii.gz$');
    
    mean_ls = spm_select('FPListRec',...
        fullfile(sub_folder), ...
        '^meanURS.*.nii$');
    
    resliced_mean_ls = spm_select('FPListRec',...
        fullfile(source_folder, ['Subject_' subs{i_sub}], 'BetaMapping'), ...
        '^r4meanURS.*.nii.gz$');
       
    rp_ls = spm_select('FPListRec',...
        fullfile(sub_folder), ...
        '^rp.*.txt$');
    
    task = '_task-audiovisualattention_run-';
    
    for i_func = 1:size(func_ls,1)
        
        if i_func<3
            ses = '1';
        else
            ses = '2';
        end
        
        run_nb = 2 - mod(i_func, 2);
        
        file = ['sub-' subs{i_sub} '_ses-' ses task sprintf('0%1.0f', run_nb) '_bold'];
        
        % transfer the func files
        SOURCE = deblank(func_ls(i_func,:));
        DESTINATION = fullfile( dest_folder, ['ses-' ses], 'func',...
            ['UR' file '.nii.gz'] );
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
        
        % transfer the rp files
        SOURCE = deblank(rp_ls(i_func,:));
        DESTINATION = fullfile( dest_folder, ['ses-' ses], 'func',...
            ['rp' file '.txt'] );
%         [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
        
    end
    
    % transfer the mean file
    
    file = ['sub-' subs{i_sub} '_ses-1' task '01_bold'];
    
    SOURCE = mean_ls;
    DESTINATION = fullfile( dest_folder, 'ses-1', 'func',...
        ['meanUR' file '.nii'] );
%     [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    
    SOURCE = resliced_mean_ls;
    DESTINATION = fullfile( dest_folder, 'ses-1', 'func',...
        ['r4meanUR' file '.nii.gz'] );
    [SUCCESS,MESSAGE,MESSAGEID] = movefile(SOURCE,DESTINATION);
    
%     copyfile(DESTINATION, fullfile(source_folder, ['Subject_' subs{i_sub}], 'BetaMapping'))
    
end
