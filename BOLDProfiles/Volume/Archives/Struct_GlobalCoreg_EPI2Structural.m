clear;  clc;

% StartFolder = fullfile(pwd, '..', '..');
StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

addpath(genpath('/data/AV_Integration_2/SubFun'))

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


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

flags = struct(...
    'sep', [4 2 1 .5] , ...
    'params',  [0 0 0  0 0 0], ...
    'cost_fun', 'nmi', ...
    'tol', [repmat(0.001, 1, 3), repmat(0.0005, 1, 3), repmat(0.005, 1, 3), repmat(0.0005, 1, 3)], ...
    'fwhm', [7,7], ...
    'graphics', ~spm('CmdLine'));


%% Define folders, number of runs, scans per run...
for SubjInd = 9:size(SubjectList,1)
    
    SubjID=SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);

    ImagesFiles2Process = {};
    TargetScan = []; %#ok<NASGU>
    SourceScan = []; %#ok<NASGU>
    
    cd(fullfile(SubjectFolder, 'Structural','CBS'))
    Struct = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii']);
    TargetScan = fullfile(SubjectFolder, 'Structural','CBS', Struct.name);
    
    
%     AnalysisFolder = fullfile(SubjectFolder, 'Transfer');
    AnalysisFolder = fullfile(SubjectFolder, 'Nifti', 'NoMoCo');
    
    cd(AnalysisFolder)
    
    MeanImg = dir('mean*.nii');
    SourceScan = fullfile(AnalysisFolder, MeanImg.name);
    
%     ImagesFilesList = dir('con*.nii');
%     
%     for FileInd = 1:length(ImagesFilesList)
%         ImagesFiles2Process{end+1,1} = fullfile(AnalysisFolder, strcat(ImagesFilesList(FileInd).name));
%     end
    
    spm_coreg_reorient_save(TargetScan,SourceScan,ImagesFiles2Process,flags);
        
    cd(StartFolder)
end


