function Analysis_MVPA_02_MVPA_Vertices
clc; clear;

StartDirectory = pwd;

addpath(fullfile(StartDirectory, 'SubFun'))


SubjectList = [...
%     '02';...
%     '03';...
    '04';...
    %'06';...
    %'07';...
%     '08';...
%     '09';...
%     '11';...
%     '12';...
%     '13';...
    %'14';...
%     '15';...
%     '16'
    ];

FoldersNames = {...
%     1:4;...
%     1:4;...
    1:4;...
    %1:2;...
    %1:4;...
%     1:4;...
%     1:4;...
%     1:4;...
%     1:4;...
%     1:4;...
    %1:2;...
%     1:4;...
%     1:4;...
    };



FFX = {'0'};


% --------------------------------------------------------- %
%                            ROIs                           %
% --------------------------------------------------------- %
ROIs_Ori = {...
    'TE1.0', [1],  {}, {}; ...
    'TE1.1', [2],  {}, {}; ...
    'TE1.2', [3],  {}, {}; ...
    'TE', [1 2 3],  {}, {}; ...
    'STG', [1 2 3 4],  {}, {}; ...
    'hA', [4],  {}, {}; ...
    
    'V1', 11,  {}, {}; ...
    'V2', 12,  {}, {}; ...
    'V3', [13],  {}, {}; ...
    'V4', [14],  {}, {}; ...
    'V5', 15,  {}, {}; ...
    
    };
NbLayers = 7;
NbLayers = NbLayers+1;

NbWorkers = 4;

% Number of beta images per condition per session (in case each block is
% modeled by a different regressor).
BlockPerCond = 3;

% Options for the SVM
opt.verbose = 0;
opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

% Saves an image with...
% ...the voxel kept for the final model
opt.output.mask = 0;
% ...the weight of each of those voxel
opt.output.weight = 0;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'}; %#ok<NASGU>


% --------------------------------------------------------- %
%              Classes and associated conditions            %
% --------------------------------------------------------- %
Class(1) = struct('name', 'A Stim', 'cond', cell(1), 'nbetas', 2*BlockPerCond);
Class(end).cond = {'A Stim - Auditory Attention' 'A Stim - Visual Attention'};

Class(2) = struct('name', 'V Stim', 'cond', cell(1), 'nbetas', 2*BlockPerCond);
Class(end).cond = {'V Stim - Auditory Attention' 'V Stim - Visual Attention'};

Class(3) = struct('name', 'AV Stim', 'cond', cell(1), 'nbetas', 2*BlockPerCond);
Class(end).cond = {'AV Stim - Auditory Attention' 'AV Stim - Visual Attention'};

Class(4) = struct('name', 'A+V Stim', 'cond', cell(1), 'nbetas', 4*BlockPerCond);
Class(end).cond = {'A Stim - Auditory Attention' 'A Stim - Visual Attention' ...
    'V Stim - Auditory Attention' 'V Stim - Visual Attention'};


Class(5) = struct('name', 'A Att', 'cond', cell(1), 'nbetas', 3*BlockPerCond);
Class(end).cond = {'A Stim - Auditory Attention' 'V Stim - Auditory Attention' 'AV Stim - Auditory Attention'};

Class(6) = struct('name', 'V Att', 'cond', cell(1), 'nbetas', 3*BlockPerCond);
Class(end).cond = {'A Stim - Visual Attention' 'V Stim - Visual Attention' 'AV Stim - Visual Attention'};



Class(7) = struct('name', 'A Att A Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'A Stim - Auditory Attention'};

Class(8) = struct('name', 'V Att A Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'A Stim - Visual Attention'};


Class(9) = struct('name', 'A Att V Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'V Stim - Auditory Attention'};

Class(10) = struct('name', 'V Att V Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'V Stim - Visual Attention'};


Class(11) = struct('name', 'A Att AV Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'AV Stim - Auditory Attention'};

Class(12) = struct('name', 'V Att AV Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
Class(end).cond = {'AV Stim - Visual Attention'};


% --------------------------------------------------------- %
%          Data pre-processing and SVM parameters           %
% --------------------------------------------------------- %
% Normalization and scaling
opt.scaling.idpdt = 0;

opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 0;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;

% Feature selection (FS)
opt.fs.threshold = 0.75;
opt.fs.type = 'ttest2';

% Recursive feature elminiation (RFE)
opt.rfe.threshold = 0.01;
opt.rfe.nreps = 20;

% SVM C/nu parameters and default arguments
opt.svm.machine = 'C-SVC';
if strcmp(opt.svm.machine, 'C-SVC')
    opt.svm.log2c = 2.^(-10:2:10);
    opt.svm.dargs = '-s 0';
elseif strcmp(opt.svm.machine, 'nu-SVC')
    opt.svm.nu = [0.0001 0.001 0.01 0.1:0.1:1];
    opt.svm.dargs = '-s 1';
end

opt.svm.kernel = 0;
if opt.svm.kernel
    % should be implemented
else
    opt.svm.dargs = [opt.svm.dargs ' -t 0 -q']; % inherent linear kernel, quiet mode
end

% Randomization options
if opt.permutation.test;
    opt.permutation.nreps = 101; % #repetitions for permutation test
else
    opt.permutation.nreps = 1;
end

% Learning curve
% #repetitions for session subsampling if needed
opt.session.subsample.nreps = 30;


MatlabVer = version('-release');
if str2double(MatlabVer(1:4))>2013
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(NbWorkers);
    end
else
    if matlabpool('size') == 0 %#ok<*DPOOL>
        matlabpool(NbWorkers)
    elseif matlabpool('size') ~= NbWorkers
        matlabpool close
        matlabpool(NbWorkers)
    end
end

for FS = 0
    
    opt.fs.do = FS;
    
    for RFE = 0
        
        opt.rfe.do = RFE;
        
        for Idpdt = 1
            
            opt.scaling.idpdt = Idpdt;
            
            for Norm = 2
                
                switch Norm
                    case 1
                        opt.scaling.img.eucledian = 1;
                        opt.scaling.img.zscore = 0;
                        opt.scaling.feat.mean = 0;
                        opt.scaling.feat.range = 0;
                        opt.scaling.feat.sessmean = 1;
                    case 2
                        opt.scaling.img.eucledian = 1;
                        opt.scaling.img.zscore = 0;
                        opt.scaling.feat.mean = 0;
                        opt.scaling.feat.range = 1;
                        opt.scaling.feat.sessmean = 0;
                    case 3
                        opt.scaling.img.eucledian = 1;
                        opt.scaling.img.zscore = 0;
                        opt.scaling.feat.mean = 1;
                        opt.scaling.feat.range = 0;
                        opt.scaling.feat.sessmean = 0;
                    case 4
                        opt.scaling.img.eucledian = 0;
                        opt.scaling.img.zscore = 1;
                        opt.scaling.feat.mean = 0;
                        opt.scaling.feat.range = 0;
                        opt.scaling.feat.sessmean = 1;
                    case 5
                        opt.scaling.img.eucledian = 0;
                        opt.scaling.img.zscore = 1;
                        opt.scaling.feat.mean = 0;
                        opt.scaling.feat.range = 1;
                        opt.scaling.feat.sessmean = 0;
                    case 6
                        opt.scaling.img.eucledian = 0;
                        opt.scaling.img.zscore = 1;
                        opt.scaling.feat.mean = 1;
                        opt.scaling.feat.range = 0;
                        opt.scaling.feat.sessmean = 0;
                end
                
                
                for iFFX = 1 %length(FFX)
                    
                    SaveSufix = '_results_surf';
                    if opt.fs.do
                        SaveSufix = [SaveSufix '_FS']; %#ok<*AGROW>
                    end
                    if opt.rfe.do
                        SaveSufix = [SaveSufix '_RFE'];
                    end
                    if opt.permutation.test
                        SaveSufix = [SaveSufix '_Perm'];
                    end
                    if opt.session.curve
                        SaveSufix = [SaveSufix '_Lear'];
                    end
                    if opt.scaling.idpdt
                        SaveSufix = [SaveSufix '_Idpdt'];
                    end
                    if opt.scaling.img.zscore
                        SaveSufix = [SaveSufix '_ZScore'];
                    end
                    if opt.scaling.img.eucledian
                        SaveSufix = [SaveSufix '_Eucl'];
                    end
                    if opt.scaling.feat.mean
                        SaveSufix = [SaveSufix '_MeanCent'];
                    end
                    if opt.scaling.feat.range
                        SaveSufix = [SaveSufix '_Range'];
                    end
                    if opt.scaling.feat.sessmean
                        SaveSufix = [SaveSufix '_SessMeanCent'];
                    end
                    
                    SaveSufix = [SaveSufix '_FWHM_' FFX{iFFX} '.mat'];
                    
                    
                    
                    for SubjInd = 1:size(SubjectList,1)
                        
                        % --------------------------------------------------------- %
                        %                        Subject data                       %
                        % --------------------------------------------------------- %
                        SubjID = SubjectList(SubjInd,:);
                        
                        fprintf('\n\nAnalysing subject %s\n', SubjID)
                        
                        SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
                        
                        NbRuns = length(FoldersNames{SubjInd});
                        
                        % If we want to have a learning curve
                        if opt.session.curve
                            % #sessions over which to run the learning curve
                            opt.session.nsamples = [6:2:NbRuns];
                        else
                            opt.session.nsamples = NbRuns;
                        end
                        
                        AnalysisFolder = fullfile(SubjectFolder, 'Transfer', 'SVM');
                        
                        if ~isdir(AnalysisFolder)
                            mkdir(AnalysisFolder)
                        end
                        
                        ROIs = ROIs_Ori;
                        
                        % --------------------------------------------------------- %
                        %                     Analysis to perform                   %
                        % --------------------------------------------------------- %
                        SVM(1) = struct('name', 'A Stim VS V Stim', 'class', [1 2]);
                        SVM(end+1) = struct('name', 'A Stim VS AV Stim', 'class', [1 3]);
                        SVM(end+1) = struct('name', 'V Stim VS AV Stim', 'class', [2 3]);
                        SVM(end+1) = struct('name', 'A Att VS V Att', 'class', [5 6]);
                        SVM(end+1) = struct('name', 'A Stim(A Att VS V Att)', 'class', [7 8]);
                        SVM(end+1) = struct('name', 'V Stim(A Att VS V Att)', 'class', [9 10]);
                        SVM(end+1) = struct('name', 'AV Stim(A Att VS V Att)', 'class', [11 12]);
                        
                        
                        
                        %% Get the vertices indices for the ROIs
                        
                        fprintf(' Reading ROI vertices indices\n')
                        
                        if exist(fullfile(AnalysisFolder, 'ROIs.mat'), 'file')
                            
                            load(fullfile(AnalysisFolder, 'ROIs.mat'));
                            
                        else
                            
                            % LEFT
                            % Manually defined A1
                            %         LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
                            %             ['A1_T1_' SubjID '*lcr*']));
                            %         [~, ~, ROI_mapping_A1_L] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
                            %             LogFileList.name), 0, 1); clear LogFileList
                            %
                            %         ROIs{1,3} = find(ROI_mapping_A1_L==ROIs{1,2}); % Stores the vertices indices
                            
                            % Cytoarchitechtonicly defined ROIs
                            LogFileList = dir(fullfile(StartDirectory, 'Subjects_Data', ['T1_' SubjID '*lcr*.vtk']));
                            [~, ~, ROI_mapping_All_L] = read_vtk(fullfile(StartDirectory, 'Subjects_Data', LogFileList.name),...
                                0, 1); clear LogFileList
                            
                            for iROI = 1:size(ROIs,1)
                                temp = [];
                                for i=1:length(ROIs{iROI,2})
                                    temp = [temp find(ROI_mapping_All_L==ROIs{iROI,2}(i))]; % Stores the vertices indices
                                end
                                ROIs{iROI,3} = temp;
                            end
                            clear temp
                            
                            %         if size(ROI_mapping_All_L,2)~=size(ROI_mapping_A1_L,2)
                            %             error('VTK files for left ROIs have differing number of vertices.')
                            %         end
                            
                            % Stores the number of vertices in that surface for future sanity check
                            NbVertices(1) = size(ROI_mapping_All_L,2);
                            
                            clear ROI_mapping_All_L ROI_mapping_A1_L
                            
                            
                            % RIGHT
                            %         LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
                            %             ['A1_T1_' SubjID '*rcr*']));
                            %
                            %         if isempty(LogFileList)
                            %             ROIs{1,4} = [];
                            %         else
                            %             [~, ~, ROI_mapping_A1_R] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
                            %                 LogFileList.name), 0, 1); clear LogFileList
                            %
                            %             ROIs{1,4} = find(ROI_mapping_A1_R==ROIs{1,2}+100);
                            %         end
                            
                            LogFileList = dir(fullfile(StartDirectory, 'Subjects_Data', ['T1_' SubjID '*rcr*.vtk']));
                            [~, ~, ROI_mapping_All_R] = read_vtk(fullfile(StartDirectory, 'Subjects_Data', LogFileList.name) ...
                                , 0, 1); clear LogFileList
                            
                            for iROI = 1:size(ROIs,1)
                                temp = [];
                                for i=1:length(ROIs{iROI,2})
                                    temp = [temp find(ROI_mapping_All_R==ROIs{iROI,2}(i))];
                                end
                                ROIs{iROI,4} = temp;
                            end
                            clear temp i
                            
                            %         if exist('ROI_mapping_A1_R', 'var')
                            %             if size(ROI_mapping_All_R,2)~=size(ROI_mapping_A1_R,2)
                            %                 error('VTK files for right ROIs have differing number of vertices.')
                            %             end
                            %         end
                            
                            NbVertices(2) = size(ROI_mapping_All_R,2);
                            clear ROI_mapping_All_R ROI_mapping_A1_R SubXP
                            
                            save(fullfile(AnalysisFolder, 'ROIs.mat'), 'NbVertices', 'ROIs');
                            
                        end
                        
                        disp(ROIs)
                        
                        
                        %% Creates a cell that lists the names of the beta images as well as their column number in the design matrix
                        load(fullfile(SubjectFolder, 'Transfer', 'SPM.mat'))
                        
                        MAX = 0;
                        for i = 1:size(SPM.Sess,2)
                            MAX = max([MAX, length(SPM.Sess(i).col)]);
                        end
                        RegNumbers = nan(size(SPM.Sess,2),MAX);
                        for i = 1:size(SPM.Sess,2)
                            RegNumbers(i,1:length(SPM.Sess(i).col)) = SPM.Sess(i).col;
                        end
                        clear i MAX
                        
                        BetaNames = char(SPM.xX.name');
                        
                        BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
                            strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
                            strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
                            strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
                            strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
                            2);
                        
                        BetaOfInterest = find(BetaOfInterest);
                        
                        BetaNames(:,1:6)=[];
                        
                        clear SPM
                        
                        
                        %% Creates a dataset that lists for each beta of interest:
                        %   - its corresponding class
                        %   - the session in which it occurs
                        CV_Mat_Orig = [zeros(NbRuns*sum([Class.nbetas]), 1) ...
                            zeros(NbRuns*sum([Class.nbetas]), 1)] ;
                        
                        % For each condition of each class we figure out what is the associated
                        % regressors and in which sessions they occur.
                        irow = 1; iFile = 1;
                        for iClass=1:numel(Class)
                            for iCond=1:numel(Class(iClass).cond)
                                
                                tmp=BetaNames(BetaOfInterest,1:length(char(Class(iClass).cond(iCond))));
                                TEMP = BetaOfInterest(strcmp(Class(iClass).cond(iCond), cellstr(tmp)));
                                
                                for i=1:length(TEMP)
                                    CV_Mat_Orig(irow,1) = iClass;
                                    [I,~] = find(TEMP(i)==RegNumbers);
                                    CV_Mat_Orig(irow,2) = I;
                                    irow = irow + 1;
                                    
                                    BetaList{iFile,1} = sprintf('%04d', TEMP(i));
                                    iFile = iFile + 1;
                                end
                                clear TEMP I i tmp
                                
                            end
                        end
                        clear irow iClass iCond BetaNames iFile RegNumbers BetaOfInterest
                        
                        
                        %% Read features
                        fprintf(' Reading features\n')
                        
                        if exist(fullfile(AnalysisFolder, 'AllFeatures.mat'), 'file')
                            
                            load(fullfile(AnalysisFolder, 'AllFeatures.mat'));
                            
                        else
                            
                            Betas = dir(fullfile(SubjectFolder, 'BetaMapping', 'exp-0000', 'Beta*lcr.vtk'));
                            Betas = char({Betas.name}');
                            Betas = Betas(:,6:9);
                            Betas = cellstr(Betas);
                            
                            tmp = dir(fullfile(SubjectFolder, 'BetaMapping', 'exp-0000', 'Beta*lcr.vtk'));
                            for i=1:numel(tmp)
                                FilesCell{i,1}=fullfile(SubjectFolder, 'BetaMapping', 'exp-0000', tmp(i).name);
                            end
                            tmp = dir(fullfile(SubjectFolder, 'BetaMapping', 'exp-0000', 'Beta*rcr.vtk'));
                            for i=1:numel(tmp)
                                FilesCell{i,2}=fullfile(SubjectFolder, 'BetaMapping', 'exp-0000', tmp(i).name);
                            end

                            ReadFeatVTK(AnalysisFolder, NbLayers, ROIs, NbVertices, BetaList, Betas, FilesCell, NbWorkers)
                            
                            load(fullfile(AnalysisFolder, 'AllFeatures.mat'));
                            
                        end
                        
                        
                        %% Run cross-validation for each model and ROI
                        for iSVM=1:numel(SVM)
                            
                            for iROI=1:size(ROIs,1)
                                
                                SVM(iSVM).ROI(iROI).name = ROIs{iROI,1};
                                SVM(iSVM).ROI(iROI).size(1,1) = numel(ROIs{iROI,3});
                                SVM(iSVM).ROI(iROI).size(1,2) = numel(ROIs{iROI,4});
                                
                                fprintf('Analysing subject %s\n', SubjID)
                                fprintf(' Running SVM:  %s\n', SVM(iSVM).name)
                                fprintf('  Running ROI:  %s\n', SVM(iSVM).ROI(iROI).name)
                                fprintf('  Number of vertices before FS/RFE: lh %i rh %i\n', SVM(iSVM).ROI(iROI).size)
                                
                                
                                for iLayer = 1:numel(NbLayers)
                                    
                                    fprintf('   Running on %i layers\n', NbLayers(iLayer))
                                    
                                    FeaturesL = FeaturesAll{iROI,1,iLayer}; %#ok<USENS>
                                    LogFeatL = ~any(isnan(FeaturesL));
                                    FeaturesLayersL = repmat(NbLayers(iLayer):-1:1, 1, size(FeaturesL,2)/NbLayers(iLayer));
                                    
                                    FeaturesR = FeaturesAll{iROI,2,iLayer};
                                    LogFeatR = ~any(isnan(FeaturesR));
                                    FeaturesLayersR = repmat(NbLayers(iLayer):-1:1, 1, size(FeaturesR,2)/NbLayers(iLayer));
                                    
                                    FeaturesBoth = [FeaturesL FeaturesR];
                                    LogFeatBoth = [LogFeatL LogFeatR];
                                    FeaturesLayersBoth = [FeaturesLayersL FeaturesLayersR];
                                    
                                    % RNG init
                                    rng('default');
                                    opt.seed = rng;
                                    
                                    Class_Acc = struct('TotAcc', []);
                                    
                                    %% Subsample sessions for the learning curve (otherwise take all of them)
                                    for NbSess2Incl = opt.session.nsamples
                                        
                                        % All possible ways of only choosing X sessions of the total
                                        CV_id = nchoosek(1:NbRuns, NbSess2Incl);
                                        CV_id = CV_id(randperm(size(CV_id, 1)),:);
                                        
                                        % Limits the number of permutation if too many
                                        if size(CV_id, 1) > opt.session.subsample.nreps
                                            CV_id = CV_id(1:opt.session.subsample.nreps,:);
                                        end
                                        
                                        %% Subsampled sessions loop
                                        for iSubSampSess=1:size(CV_id, 1)
                                            
                                            % Permutation test
                                            if NbSess2Incl < NbRuns
                                                NbPerm = 1;
                                                fprintf('\n    Running learning curve with %i sessions\n', NbSess2Incl)
                                                fprintf('     %i of %i \n\n', iSubSampSess, size(CV_id, 1))
                                            else
                                                fprintf('    Running analysis with all sessions\n\n')
                                                NbPerm = opt.permutation.nreps;
                                            end
                                            
                                            %%
                                            for iPerm=1:NbPerm
                                                
                                                CV_Mat = CV_Mat_Orig;
                                                
                                                %% Permute class within sessions when all sessions are included
                                                if iPerm > 1
                                                    for iRun=1:max(CV_Mat.session)
                                                        temp = CV_Mat((all([ismember(CV_Mat(:,1), SVM(iSVM).class), ...
                                                            ismember(CV_Mat(:,2), iRun)], 2)),1);
                                                        
                                                        CV_Mat(:,1) = CV_Mat(temp(randperm(length(temp))),1);
                                                    end
                                                end
                                                clear temp
                                                
                                                %% Leave-one-run-out (LORO) cross-validation
                                                parfor iCV=1:size(CV_id, 2)
                                                    
                                                    fprintf('     Cross-Validation %i of %i', iCV, size(CV_id, 2))
                                                    
                                                    TestSess = []; %#ok<NASGU>
                                                    TrainSess = []; %#ok<NASGU>
                                                    
                                                    % Separate training and test sessions
                                                    [TestSess, TrainSess] = deal(false(size(1:NbRuns)));
                                                    TestSess(CV_id(iSubSampSess,iCV)) = 1;
                                                    TrainSess(setdiff(CV_id(iSubSampSess,:), CV_id(iSubSampSess,iCV))) = 1;
                                                    
                                                    [acc_layer, results_layer, results, ~] = RunSVM(SVM, FeaturesBoth, LogFeatBoth, FeaturesLayersBoth, CV_Mat, TrainSess, TestSess, opt, iSVM);
                                                    
                                                    TEMP(iCV,1, iLayer).layers.results = {results_layer};
                                                    TEMP(iCV,1, iLayer).layers.acc = acc_layer;
                                                    TEMP(iCV,1,iLayer).results = {results};
                                                    TEMP(iCV,1,iLayer).acc = mean(results.pred==results.label);
                                                    
                                                    fprintf('\n')
                                                    
                                                end
                                                
                                                for iCV=1:size(CV_id, 2)
                                                    for iSide=1:size(TEMP,2)
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,iSide, iLayer).layers.results =...
                                                            TEMP(iCV,iSide, iLayer).layers.results;
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,iSide, iLayer).layers.acc = ...
                                                            TEMP(iCV,iSide, iLayer).layers.acc;
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,iSide,iLayer).results = ...
                                                            TEMP(iCV,iSide,iLayer).results;
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,iSide,iLayer).acc = ...
                                                            TEMP(iCV,iSide,iLayer).acc;
                                                    end
                                                end
                                                
                                            end % iPerm=1:NbPerm
                                            %clear iPerm
                                            
                                        end % iSubSampSess=1:size(CV_id, 1)
                                        %clear iSubSampSess
                                        
                                    end % NbSess2Incl = opt.session.nsamples
                                    %clear NbSess2Incl
                                    
                                end % iLayer = 1:numel(NbLayers)
                                clear iLayer
                                
                                %% Calculate prediction accuracies
                                for iLayer = 1:numel(NbLayers)
                                    for Side = 1
                                        Class_Acc.TotAcc(Side,iLayer) = ...
                                            nanmean([SVM(iSVM).ROI(iROI).session(end).rand.perm.CV(:,Side,iLayer).acc]);
                                        for iCV=1:size(CV_id, 2)
                                            temp(:,:,iCV) = SVM(iSVM).ROI(iROI).session(end).rand.perm.CV(iCV,Side,iLayer).layers.acc;
                                        end
                                        Class_Acc.TotAccLayers{Side,iLayer} = nanmean(temp,3);
                                        temp = [];
                                    end
                                    
                                    % Display some results
                                    fprintf('\n   Accuracy %i layers\n', NbLayers(iLayer));
                                    disp(Class_Acc.TotAcc(:,iLayer))
                                    disp(Class_Acc.TotAccLayers{1,iLayer})
                                    
                                end
                                
                                % Save data
                                Results = SVM(iSVM).ROI(iROI);
                                SaveResults(AnalysisFolder, Results, opt, SVM, iSVM, iROI, SaveSufix)
                                
                            end % iSubROI=1:numel(SVM(iSVM).ROI)
                            clear Mask Features
                            
                        end % iSVM=1:numel(SVM)
                        clear iSVM SVM FeaturesAll
                        
                    end
                end
            end
        end
    end
end

if ~isempty(gcp)
    delete(gcp)
end

end

function SaveResults(SaveDir, Results, opt, SVM, iSVM, iROI, SaveSufix) %#ok<INUSL>

save(fullfile(SaveDir, ['SVM_' SVM(iSVM).name '_ROI_' SVM(iSVM).ROI(iROI).name SaveSufix]), 'Results', 'opt', '-v7.3');

end

function [acc_layer, results_layer, results, weight] = RunSVM(SVM, Features, LogFeat, FeaturesLayers, CV_Mat, TrainSess, TestSess, opt, iSVM)

if isempty(Features) || all(Features(:)==Inf)
    
    warning('Empty ROI')
    
    acc_layer = NaN;
    results = struct();
    results_layer = struct();
    weight = [];
    
else
    
    [acc_layer, weight, results_layer] = machine_SVC_layers(SVM(iSVM), ...
        Features(:,LogFeat), FeaturesLayers(:,LogFeat), CV_Mat, TrainSess, TestSess, opt);
    
    if opt.verbose
                fprintf('\n       Running on all layers.')
    end
    
    results = machine_SVC(SVM(iSVM), Features(:,LogFeat), CV_Mat, TrainSess, TestSess, opt);
    
    %     weight(1,:,2) = [results.w / sum(results.w)]'; %#ok<NBRAK>
    %
    %     results = rmfield(results, 'w');
    
    
end

end
