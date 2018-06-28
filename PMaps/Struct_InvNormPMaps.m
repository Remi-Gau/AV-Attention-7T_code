clear; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>

% Folders definitions
RootFolder = fullfile(pwd,'..','..');

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


for SubjInd = 9 %size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(RootFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    SegmantationFolder = fullfile(SubjectFolder, 'Structural', 'SPM_Segmentation');
    
    PMap_Folder = fullfile(SubjectFolder, 'PMap');
    mkdir(PMap_Folder)
    
    DefField = fullfile(SegmantationFolder, 'iy_UNI.nii');
    
    %     cd(fullfile(RootFolder, 'PMap', 'Cyt'))
    %     copyfile('m*.nii', PMap_Folder)
    
%     cd(fullfile(RootFolder, 'PMap', 'Ret'))
%     copyfile('perc*.nii', PMap_Folder)
    
    cd(fullfile(RootFolder, 'PMap', 'BT'))
    copyfile('1.25*.nii', PMap_Folder)
    
    cd(PMap_Folder)
    
    ImagesFiles2Process={};
    %     delete('rrw*.nii')
    %     delete('w*.nii')
    ImgLs = dir('1.25*.nii');
    
%     tmp =  dir('perc*.nii');
%     ImgLs = cat(1,ImgLs,tmp);
    
    clear tmp
    for ImageInd=1:length(ImgLs)
        ImagesFiles2Process{ImageInd,1} = fullfile(PMap_Folder, strcat(ImgLs(ImageInd).name, ',1')); %#ok<*SAGROW>
    end
    
    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    DefField
    
    ImagesFiles2Process
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {DefField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ImagesFiles2Process;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN*ones(2,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = NaN*ones(1,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    save (fullfile(SubjectFolder, 'InvNormPMap_matlabbatch.mat'), 'matlabbatch');
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
end

