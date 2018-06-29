clear; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

% Folders definitions
RootFolder = fullfile(pwd, '..', '..');
addpath(genpath(fullfile(RootFolder, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
%     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
%     '14';...
    '15';...
    '16'
    ];


for SubjInd = 1:size(SubjectList,1)
    
    clear matlabbatch
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>

    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    SegmantationFolder = fullfile(SubjectFolder, 'Structural', 'SPM_Segmentation');
    
    ROI_Folder = fullfile(SubjectFolder, 'ROI_MNI');
    mkdir(ROI_Folder)
        
%     mkdir(SegmantationFolder)   
%     copyfile(fullfile('/','media','rxg243','BackUp2','AV_Integration_7T_2', ...
%         'Subjects_Data', ['Subject_' SubjID], 'FFX', 'Structural','*.nii'), SegmantationFolder)
%     
%     copyfile(fullfile('/','media','rxg243','BackUp2','AV_Integration_7T_2', ...
%         'Subjects_Data', ['Subject_' SubjID], 'FFX', 'Structural','*.mat'), SegmantationFolder)
    
    DefField = fullfile(SegmantationFolder, 'iy_UNI.nii');
        
    cd(fullfile(RootFolder, 'ROI_MNI'))

%     copyfile('ProbRet*v.nii', ROI_Folder)
%     copyfile('ProbRet*d.nii', ROI_Folder)
%     copyfile('*S1*.*', ROI_Folder)
%     copyfile('TE_1.*_p>.6_*.*', ROI_Folder)
    
    copyfile('*BrainTome.nii', ROI_Folder)
    
    
    cd(ROI_Folder)
    ImagesFiles2Process={};
    
%     filesdir = dir('*S1*.*');
%     tmp = dir('ProbRet*d.nii');
%     tmp2 = dir('ProbRet*v.nii');

%      filesdir = dir('TE_1.*_p>.6_*.*');
     
     filesdir = dir('*BrainTome.nii');
    
%     filesdir = [filesdir;tmp;tmp2]; %#ok<AGROW>
    
    clear tmp
    for ImageInd=1:length(filesdir)
        ImagesFiles2Process{ImageInd,1} = fullfile(ROI_Folder, strcat(filesdir(ImageInd).name, ',1')); %#ok<*SAGROW>
    end

    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {DefField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ImagesFiles2Process;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN*ones(2,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = NaN*ones(1,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
        
    save (fullfile(SubjectFolder, 'InvNorm_matlabbatch.mat'), 'matlabbatch');
    
    spm_jobman('run',matlabbatch)   
    
    cd (RootFolder)
    
end

