clear; clc

% Folders definitions
StartDirectory = fullfile(pwd, '..', '..');

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


ROIs= {...
    'A1_surf'
    'PT_surf_thres'
    'V1_surf_thres'
    'V2-3_surf_thres'
    };


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    ROI_Folder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    ActMaskFolder = fullfile('/media/rxg243/BackUp2/AV_Integration_7T_2', 'Subjects_Data', ['Subject_' SubjID], 'FFX_Block_Smooth'); 
    
    cd(ActMaskFolder)
    ActMaskHdr = dir('rspmT_*Baseline_thres_0.1.nii');
    ActMaskHdr = spm_vol(char({ActMaskHdr.name}'));
    ActMaskVol = spm_read_vols(ActMaskHdr);
    
    cd(ROI_Folder)
    for iROI=1:numel(ROIs)
        ROI = spm_read_vols(spm_vol([ROIs{iROI} '.nii']));
  
        for iMask=1:numel(ActMaskHdr)
            Vol = all(cat(4,ROI,ActMaskVol(:,:,:,iMask)),4);
            Hdr = ActMaskHdr(1);
            Hdr.fname = [ROIs{iROI} '_' strrep(ActMaskHdr(iMask).fname(12:end-4), ' ', '_') '.nii'];  
            spm_write_vol(Hdr, Vol)
            if sum(Vol(:))==0
                warning('Intersection of %s and %s : No voxel left', ROIs{iROI}, ActMaskHdr(iMask).fname)
            end
        end
    end

    cd (StartDirectory)
    
end

