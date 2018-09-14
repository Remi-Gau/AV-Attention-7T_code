% A small script that tries to estiame the difference of doing the realign
% with and witouth the field map. Also tries to estimate the change in
% smoothness between the 2.

clear; clc;

% StartDirectory = fullfile(pwd, '..','..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';
cd (StartDirectory)

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
    
    AnalysisFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID] ...
        , 'Nifti', 'NoMoCo', '01');
    
    cd(AnalysisFolder)
    
    RealAndUnwarpWithFM = spm_read_vols(spm_vol(spm_select('FPList',pwd,'^meanURS.*\.nii$')));
    RealAndUnWarpNoFM = spm_read_vols(spm_vol(spm_select('FPList',pwd,'^meanuS.*\.nii$')));
    RealThenUnWarp = spm_read_vols(spm_vol(spm_select('FPList',pwd,'^meanu2S.*\.nii$')));
    
    
    hdr = spm_vol(spm_select('FPList',pwd,'^meanURS.*\.nii$'));
    hdr.fname = fullfile(pwd,'diff_mean_1.nii');
    spm_write_vol(hdr,  (RealAndUnwarpWithFM-RealThenUnWarp).^2)
    
    hdr.fname = fullfile(pwd,'diff_mean_2.nii');
    spm_write_vol(hdr,  (RealAndUnwarpWithFM-RealAndUnWarpNoFM).^2)
    
    
    hdr = spm_vol(spm_select('FPList',pwd,'^meanURS.*\.nii$'));
    mask = RealAndUnwarpWithFM>spm_global(hdr);
    hdr_mask = hdr;
    hdr_mask.fname = fullfile(pwd,'mask_tmp.nii');
    spm_write_vol(hdr_mask, mask)
    Smoothness(1:3,1,SubjInd) = spm_est_smoothness(spm_select('FPList',pwd,'^meanURS.*\.nii$'),spm_select('FPList',pwd,'^mask.*\.nii$'));  %#ok<*SAGROW>
    
    hdr = spm_vol(spm_select('FPList',pwd,'^meanu2S.*\.nii$'));
    mask = RealThenUnWarp>spm_global(hdr);
    spm_write_vol(hdr_mask, mask)
    Smoothness(1:3,2,SubjInd) = spm_est_smoothness(spm_select('FPList',pwd,'^meanu2S.*\.nii$'),spm_select('FPList',pwd,'^mask.*\.nii$'));
    
    hdr = spm_vol(spm_select('FPList',pwd,'^meanuS.*\.nii$'));
    mask = RealAndUnWarpNoFM>spm_global(hdr);
    spm_write_vol(hdr_mask, mask)
    Smoothness(1:3,3,SubjInd) = spm_est_smoothness(spm_select('FPList',pwd,'^meanuS.*\.nii$'),spm_select('FPList',pwd,'^mask.*\.nii$'));
    
end