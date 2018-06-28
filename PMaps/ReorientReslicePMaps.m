clear; clc;

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


StartDirectory = fullfile(pwd,'..','..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))


SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';
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


for SubjInd =  9 %1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
                            
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);                    
    
    AnalysisFolder = fullfile(SubjectFolder, 'PMap');

    %% Reorient ROIs
    cd(AnalysisFolder)
    
    Files2Reorient = {};
    
    ImgLs = dir('w1.25*.nii');  
%     ImgLs = dir('wperc*.nii');  

    for iROI=1:numel(ImgLs)
        Files2Reorient{end+1,1} = fullfile(AnalysisFolder,ImgLs(iROI).name); %#ok<*SAGROW>
    end
    clear tmp i
    
%     cd(fullfile(StartDirectory, 'Subjects_Data', 'Transfer'));
    cd(fullfile('/media/rxg243/BackUp2/AV_Integration_7T_2','Subjects_Data', ['Subject_' SubjID], 'Transfer'));
    
    MatFiles = dir('*mat');
    MatFiles = char(MatFiles.name);
    
%     Reorient = strcmp(cellstr(MatFiles(:,1:14)),'ReorientMatrix');
%     if any(Reorient)
%         if sum(Reorient)>1
%             warning('Images have been reoriented more than once.')
%         else
%             File2Load = find(Reorient);
%             for iFile=1:length(File2Load)
%                 load(deblank(MatFiles(File2Load(iFile),:)), 'M')
%                 MatFiles(File2Load(iFile),:)
%                 M
% %                 spm_reorient(Files2Reorient, M)
%             end
%         end
%     end
%     clear Reorient M
%     
%     Coregister = strcmp(cellstr(MatFiles(:,1:11)),'CoregMatrix');
%     if any(Coregister)
%         if sum(Coregister)>1
%             warning('Images have been coregistered more than once.')
%         else
%             File2Load = find(Coregister);
%             for iFile=1:length(File2Load)
%                 load(deblank(MatFiles(File2Load(iFile),:)), 'M')
%                 MatFiles(File2Load(iFile),:)
%                 M
% %                 spm_reorient(Files2Reorient, M)
%             end
%         end
%     end
    
    clear Files2Reorient M
    
    %% Reslice
    cd(AnalysisFolder)

    Files2Reslice = {};
    matlabbatch = {};
    
    ROIFiles = dir('w1.25*.nii'); 
%         ROIFiles = dir('wperc*.nii');  

    for i=1:numel(ROIFiles)
        Files2Reslice{end+1,1} = fullfile(AnalysisFolder, [ROIFiles(i).name ',1']);
    end

    matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID], 'Structural', ...
        'CBS', ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound.nii'])};
    
    Files2Reslice
    
    matlabbatch{1}.spm.spatial.coreg.write.source = Files2Reslice;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'rr';

    save (strcat('ResliceROIs_Subject_', SubjID, '_jobs.mat'));
    
    spm_jobman('run', matlabbatch)


end