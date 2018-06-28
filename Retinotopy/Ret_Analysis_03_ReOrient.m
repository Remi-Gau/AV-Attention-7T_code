%%
clear all; clc

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
    '16';...
    ];

flagsRealign = struct(...
    'sep', [1.5 1], ...
    'params',  [0 0 0  0 0 0], ...
    'cost_fun', 'nmi', ...
    'tol', [repmat(0.001, 1, 3), repmat(0.0005, 1, 3), repmat(0.005, 1, 3), repmat(0.0005, 1, 3)], ...
    'fwhm', [7,7], ...
    'graphics', ~spm('CmdLine'));

flagsReslice = struct(...
    'mask', 0, ...
    'mean', 0, ...
    'interp', 1, ...
    'which', 1, ...
    'wrap', [0 0 0], ...
    'prefix', 'r' ...
    );


% To move EPI to roughly the same space as the structural
P_translation = zeros(12,1);

P_translation(1) = 36;
P_translation(2) = -33;
P_translation(3) = -167;

TranslationMat = spm_matrix(P_translation, 'T');


%  Root folder definition
StartFolder = pwd;

cd(StartFolder)

for SubjInd = 1:size(SubjectList,1)
    
    %% Subject's Identity and folders
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    RetFolder = fullfile(SubjectFolder, 'Retinotopy');
    
    [~,~,~]=mkdir(fullfile(RetFolder, 'MIPAV'));
    
    copyfile(fullfile(RetFolder, 'FFX', 'Polar*.nii'), fullfile(RetFolder, 'MIPAV'));
    copyfile(fullfile(RetFolder, 'Nifti', '01', 'mean*.nii'), fullfile(RetFolder, 'MIPAV'));
    copyfile(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_*.nii'), fullfile(RetFolder, 'MIPAV'));
    
    
    %%
    cd(fullfile(RetFolder, 'MIPAV'));
    tmp = dir('mean*.nii');
    EPI_Hdr = spm_vol(tmp.name); 
    EPI_Vol = spm_read_vols(EPI_Hdr);

    EPI_Hdr.mat = TranslationMat * EPI_Hdr.mat;
    
    spm_write_vol(EPI_Hdr, EPI_Vol);
    
    clear tmp EPI_Hdr EPI_Vol
    
    tmp = dir('Polar*.nii');
    for i=1:length(tmp)
        EPI_Hdr = spm_vol(tmp(i).name);
        EPI_Vol = spm_read_vols(EPI_Hdr);
        
        EPI_Hdr.mat = TranslationMat * EPI_Hdr.mat;
        
        spm_write_vol(EPI_Hdr, EPI_Vol);
    end
    
    clear tmp EPI_Hdr EPI_Vol
    
    
    %%
    ImagesFiles2Process = {};
    TargetScan = []; %#ok<NASGU>
    SourceScan = []; %#ok<NASGU>
    
    tmp = dir('T1*.nii');
    TargetScan = fullfile(pwd, tmp.name);
    tmp = dir('mean*.nii');
    SourceScan = fullfile(pwd, tmp.name);
    
    ImagesFilesList = dir('Polar*.nii');
    
    for FileInd = 1:length(ImagesFilesList)
        ImagesFiles2Process{end+1,1} = fullfile(pwd, strcat(ImagesFilesList(FileInd).name));
    end
    
    spm_coreg_reorient_save(TargetScan,SourceScan,ImagesFiles2Process,flagsRealign);
    
    
    %%
    ImagesFiles2Process = {};
    
    tmp = dir('T1*.nii');
    ImagesFiles2Process = {fullfile(pwd, tmp.name);}
    tmp = dir('mean*.nii');
    ImagesFiles2Process{2,1} = fullfile(pwd, tmp.name);
    
    ImagesFilesList = dir('Polar*.nii');
    
    for FileInd = 1:length(ImagesFilesList)
        ImagesFiles2Process{end+1,1} = fullfile(pwd, strcat(ImagesFilesList(FileInd).name));
    end
    
     spm_reslice(ImagesFiles2Process,flagsReslice)
    
end