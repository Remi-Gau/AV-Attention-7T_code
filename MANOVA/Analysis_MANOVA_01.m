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

ROI(1) = struct('name', 'STGpost', 'fname', 'rrwSTG_Post_AAL.nii');
ROI(end+1) = struct('name', 'TE', 'fname', 'rrwTE_MNI.nii');
ROI(end+1) = struct('name', 'V1', 'fname', 'rrwProbRet_V1.nii');
ROI(end+1) = struct('name', 'V2-3', 'fname', 'rrwProbRet_V2-3.nii');

permute = 0;
lambda = 0;
slRadius = 5;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

for SubjInd = 1 %1:size(SubjectList,1)
    
    Analysis(1) = struct('name', 'A Stim VS V Stim', 'class1', [1 4], 'class2', [2 5]);
    Analysis(end+1) = struct('name', 'A Stim VS AV Stim', 'class1', [1 4], 'class2', [3 6]);
    Analysis(end+1) = struct('name', 'V Stim VS AV Stim', 'class1', [2 5], 'class2', [3 6]);
    Analysis(end+1) = struct('name', 'A Att VS V Att', 'class1', 1:3, 'class2', 4:6);
    
    
    %% Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\nRunning subject %s.\n\n', SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block_UpSamp');
    LayerFolder = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering');
    ROIFolder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    
    
    %% Define contrasts
    load(fullfile(AnalysisFolder, 'SPM.mat'))
    
    % Regressor numbers
    MAX = 0;
    for i = 1:size(SPM.Sess,2)
        MAX = max([MAX, length(SPM.Sess(i).col)]);
    end
    RegNumbers = nan(size(SPM.Sess,2),MAX);
    for i = 1:size(SPM.Sess,2)
        RegNumbers(i,1:length(SPM.Sess(i).col)) = SPM.Sess(i).col;
    end
    clear i MAX
    
    % Regressor names
    TEMP = char(SPM.xX.name');
    TEMP(:,1:6)=[];
    for i=1:size(TEMP,1)
        if TEMP(i,1) == ' '
            RegNames(i,:) = [TEMP(i,2:end) ' ']; %#ok<*SAGROW>
        else
            RegNames(i,:) = TEMP(i,:);
        end
    end
    clear TEMP
    
    % Create contrasts
    for iAnalysis=1:numel(Analysis)
        All = zeros(size(RegNames,1),1);
        for j=1:numel(Analysis(iAnalysis).class1)
            for k=1:3
                TEMP = strcmp([Conditions_Names{Analysis(iAnalysis).class1(j)} ' - Block ' num2str(k) '*bf(1)'], ...
                    cellstr(RegNames));
                All = All+TEMP;
            end
        end
        for j=1:numel(Analysis(iAnalysis).class2)
            for k=1:3
                TEMP = strcmp([Conditions_Names{Analysis(iAnalysis).class2(j)} ' - Block ' num2str(k) '*bf(1)'], ...
                    cellstr(RegNames));
                All = All-TEMP;
            end
        end
        Contrasts(iAnalysis,:) = All;
    end

    
    %% Get Layer labels
    LayerVol = spm_read_vols(spm_vol(fullfile(LayerFolder,['T1_' sprintf('%02.0f',NbLayers) '_Layers.nii'])));
    
    
    %% Get ROI
    cd(ROIFolder)
    ROI_Vol = spm_read_vols(spm_vol(char({ROI.fname}')));
    
    
    %% Run MANOVA
    cd(AnalysisFolder)
    for iAnalysis=1:numel(Analysis)
        Cs{iAnalysis} = Contrasts(iAnalysis,RegNumbers(1,~isnan(RegNumbers(1,:))));
    end

    s = warning('off', 'MATLAB:nearlySingularMatrix');

    for iROI = 1:numel(ROI)
        
        for iLayer = 1:NbLayers
            ROI_Layer_Vol = all( cat(4, LayerVol==iLayer, ROI_Vol(:,:,:,iROI) ), 4);
            
            cvManovaSearchlight(AnalysisFolder, slRadius, Cs, ROI_Layer_Vol, permute, lambda)
            
            % Rename files so they are not overwritten
            FileLs = dir('spmDs_C000*_P0001.nii');
            for iFile=1:numel(FileLs)
                movefile(FileLs(iFile).name, ...
                    [FileLs(iFile).name(1:end-4) '_L' num2str(iLayer) '.nii'])
            end
            
            movefile('VPSL.nii', ['VPSL_' num2str(iLayer) '.nii'])
            
        end
        
        % Puts all the information WRT VPSL in one image per ROI
        FileLs = dir('VPSL_*.nii');
        Hdr = spm_vol(char({FileLs.name}'));
        Vol = spm_read_vols(Hdr);
        Vol = nansum(Vol,4);
        Hdr = Hdr(1);
        Hdr.fname = [ROI(iROI).name '_VPSL.nii'];
        spm_write_vol(Hdr,Vol)
        
        
        % Puts all the information WRT one contrast in one image per ROI
        

    end

    warning(s.state, 'MATLAB:nearlySingularMatrix')
    
end