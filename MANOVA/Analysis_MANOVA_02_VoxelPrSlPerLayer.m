clc; clear;

StartFolder = fullfile(pwd,'..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')));

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

NbLayers = 6;
SlRadius = 8;

ROI(1) = struct('name', 'STGpost', 'fname', 'rrwSTG_Post_AAL.nii', 'HistData', cell(1));
ROI(end+1) = struct('name', 'TE', 'fname', 'rrwTE_MNI.nii', 'HistData', cell(1));
ROI(end+1) = struct('name', 'V1', 'fname', 'rrwProbRet_V1.nii', 'HistData', cell(1));
ROI(end+1) = struct('name', 'V2-3', 'fname', 'rrwProbRet_V2-3.nii', 'HistData', cell(1));

for SubjInd = 1 %1:size(SubjectList,1)
    
   
    %% Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\nRunning subject %s.\n\n', SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block_UpSamp', ['SL' num2str(SlRadius)]);
    LayerFolder = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering');
    
    %% Get Layer labels
    LayerVol = spm_read_vols(spm_vol(fullfile(LayerFolder,['T1_' sprintf('%02.0f',NbLayers) '_Layers.nii'])));
    
    
    %% Run MANOVA
    cd(AnalysisFolder)
    
    HistData = cell(NbLayers,numel(ROI));

    for iROI = 1:numel(ROI)
        
        ROI_Vol = spm_read_vols(spm_vol([ROI(iROI).name '_VPSL.nii']));

        for iLayer = 1:NbLayers

            ROI(iROI).HistData{iLayer} = ROI_Vol(LayerVol==iLayer);
                      
        end

    end
    
end

%%
figure('name', 'Voxels in searchlight per layer')
for iROI = 1:numel(ROI)
    subplot(1, numel(ROI), iROI)
    distributionPlot(ROI(iROI).HistData)
    title(ROI(iROI).name)
end