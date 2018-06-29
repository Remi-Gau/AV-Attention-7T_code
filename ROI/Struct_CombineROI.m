clear; clc

% Folders definitions
RootFolder = fullfile(pwd, '..', '..');

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

Mask.ROI(1) = struct('name', 'V1', 'fname', 'rrwProbRet_V1.nii');
Mask.ROI(end+1) = struct('name', 'V2', 'fname', 'rrwProbRet_V2.nii');
Mask.ROI(end+1) = struct('name', 'V3', 'fname', 'rrwProbRet_V3.nii');
Mask.ROI(end+1) = struct('name', 'V1v', 'fname', 'rrwProbRet_V1v.nii');
Mask.ROI(end+1) = struct('name', 'V2v', 'fname', 'rrwProbRet_V2v.nii');
Mask.ROI(end+1) = struct('name', 'V3v', 'fname', 'rrwProbRet_V3v.nii');
Mask.ROI(end+1) = struct('name', 'V1d', 'fname', 'rrwProbRet_V1d.nii');
Mask.ROI(end+1) = struct('name', 'V2d', 'fname', 'rrwProbRet_V2d.nii');
Mask.ROI(end+1) = struct('name', 'V3d', 'fname', 'rrwProbRet_V3d.nii');

Vol2Merge = {...
    'rrwProbRet_V2-3.nii' , [2:3];...
    'rrwProbRet_V1-2.nii' , [1:2];...
    'rrwProbRet_V1-2-3.nii' , [1:3];...
    'rrwProbRet_V1-2-3d.nii' , [7:9];...
    'rrwProbRet_V1-2-3v.nii' , [4:6]};

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    ROI_Folder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    
    cd(ROI_Folder)
    
    for i = 1%:size(Vol2Merge,1)
        hdr = spm_vol(char({Mask.ROI(Vol2Merge{i,2}).fname}'));
        vol = spm_read_vols(hdr);
        
        vol = any(vol,4);
        hdr = hdr(1);
        hdr.fname = Vol2Merge{i,1};
        spm_write_vol(hdr, vol)
    end
 
    
    cd (RootFolder)
    
end

