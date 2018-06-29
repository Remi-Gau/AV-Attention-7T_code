clc; clear;

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

ROI = {'A1_V_deact';'PT_V_deact';'V1_A_deact';'V2-3_A_deact'}


for SubjInd = 8:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n Subject %s\n', SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'ROI_MIPAV', 'ActDeact');
    
    for iROI=1:numel(ROI)
        
        Left = fullfile(SubjectFolder, ['exp-000' num2str(iROI-1)], ...
            ['exp-000' num2str(iROI-1) '-AA'],'MeshDataToVolume');
        
        File = spm_select('FPList', Left, ['^T1_' SubjID  '.*lcr.*_data.nii.gz$']);
%         gunzip(File)
        FileL = spm_select('FPList', Left, ['^T1_' SubjID  '.*lcr.*_data.nii$']);
        
        Right = fullfile(SubjectFolder, ['exp-000' num2str(iROI-1)], ...
            ['exp-000' num2str(iROI-1) '-BA'],'MeshDataToVolume');
        
        File = spm_select('FPList', Right, ['^T1_' SubjID  '.*rcr.*_data.nii.gz$']);
%         gunzip(File)
        FileR = spm_select('FPList', Right, ['^T1_' SubjID  '.*rcr.*_data.nii$']);
        
        
        hdr = spm_vol(cat(1,FileL,FileR));
        vol = spm_read_vols(hdr);
        
        vol = sum(vol<0,4);
        
        hdr = hdr(1);
        hdr.fname = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Transfer', 'ROI', [ROI{iROI} '.nii']);
        spm_write_vol(hdr, vol);
        
        clear FileLFileR hdr vol
        
    end
    
    
end