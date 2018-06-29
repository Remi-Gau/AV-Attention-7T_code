clear; clc;

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


StartDirectory = fullfile(pwd,'..','..','..');

addpath(fullfile(StartDirectory, 'SubFun'))


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

flags = struct(...
    'mask', 0, ...
    'mean', 0, ...
    'interp', 0, ...
    'which', 1, ...
    'wrap', [0 0 0], ...
    'prefix', 'rr' ...
    );


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
                            
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);                    

    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    load(fullfile(GLM_Folder, 'SPM.mat'));
    
    for iSess=1:numel(SPM.xX.iB)
        copyfile(fullfile(SubjectFolder, 'FFX_Block', ['beta_' sprintf('%04.0f', SPM.xX.iB(iSess)) '.nii']), ...
            fullfile(SubjectFolder, ['SessMean_' num2str(iSess) '.nii']))
    end
    
    
    

    %% Reorient ROIs
    cd(SubjectFolder)
    
    Files2Reorient = {};

    filesdir = dir('SessMean*.nii');  

    for iImg=1:numel(filesdir)
        Files2Reorient{end+1,1} = fullfile(SubjectFolder,filesdir(iImg).name); %#ok<*SAGROW>
    end
    
    cd(fullfile(SubjectFolder, 'Transfer'));
    
    MatFiles = dir('*mat');
    MatFiles = char(MatFiles.name);
    
    Reorient = strcmp(cellstr(MatFiles(:,1:14)),'ReorientMatrix');
    if any(Reorient)
        if sum(Reorient)>1
            error('Images have been reoriented more than once.')
            return
        else
            load(deblank(MatFiles(Reorient,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    clear Reorient M
    
    Coregister = strcmp(cellstr(MatFiles(:,1:11)),'CoregMatrix');
    if any(Coregister)
        if sum(Coregister)>1
            error('Images have been coregistered more than once.')
            return
        else
            load(deblank(MatFiles(Coregister,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    
    clear Files2Reorient M
    
    %% Reslice
    cd(SubjectFolder)

    Files2Reslice = {};
    matlabbatch = {};
    
    filesdir = dir('SessMean*.nii');  
    
    for i=1:numel(filesdir)
        Files2Reslice{end+1,1} = fullfile(SubjectFolder, [filesdir(i).name ',1']);
    end

    matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'rr';

    save (strcat('ResliceSessMeans_Subject_', SubjID, '_jobs.mat'));
    
    spm_jobman('run', matlabbatch)


end