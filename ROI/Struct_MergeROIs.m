clc; clear all;

StartDirectory = pwd;

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

ToMerge(1).ROI = struct('name', 'TE', 'dir', 'A', 'subroi', {{'TE_1.0'; 'TE_1.1'; 'TE_1.2'}});
ToMerge(2).ROI = struct('name', 'V3_4_d', 'dir', 'V', 'subroi', {{'V3d'; 'V4d'}});
ToMerge(3).ROI = struct('name', 'V3_4_v', 'dir', 'V', 'subroi', {{'V3v'; 'V4v'}});

HS_Sufix = {...
    '_lh','_L_MNI.nii';
    '_rh','_R_MNI.nii'};

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n Subject %s\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    for hs = 1:2
        
        for iROI=1:size(ToMerge,2)
            
            RoiFolder = fullfile(SubjectFolder, 'Analysis', 'ROI', ...
                [ToMerge(iROI).ROI.dir HS_Sufix{hs,1}]);
            
            Vol = [repmat([RoiFolder filesep 'rsub_'], size(ToMerge(iROI).ROI.subroi), 1) ...
                char(ToMerge(iROI).ROI.subroi) ...
                repmat(HS_Sufix{hs,2}, size(ToMerge(iROI).ROI.subroi), 1)];
            
            hdr = spm_vol(Vol(1,:));
            
            IMG = logical(spm_read_vols(spm_vol(Vol)));
            IMG = any(IMG,4);
            
            hdr.fname = fullfile(RoiFolder, ['rsub_' ToMerge(iROI).ROI.name HS_Sufix{hs,2}]);
            
            spm_write_vol(hdr, IMG);
            
            clear IMG hdr Vol RoiFolder
            
        end
        
    end
    
end