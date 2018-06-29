function Vertex_Analysis_06_MVPA_Layer

% Finish Subject 03 (? - ?)
% Finish Subject 11 (V Stim VS AV Stim - V2)
% Finish Subject 12 (A Stim VS V Stim - V1)
% Finish Subject 15 (A Stim VS AV Stim - V3d)

clc; clear;

StartDirectory = pwd;

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %     '06';...
    '07';...
    '08';...
    '09';...%6
    '11';...
    '12';...
    '13';...
    %     '14';...
    '15';...
    '16'
    ];

XP = [...
    'exp-0000';...
    'exp-0001';...
    'exp-0002';...
    %     '06';...
    'exp-0004';...
    'exp-0005';...
    'exp-0006';...
    'exp-0007';...
    'exp-0008';...
    'exp-0009';...
    %     '14';...
    'exp-0011';...
    'exp-0012'
    ];

FoldersNames = {...
    1:4;...
    1:4;...
    1:4;...
    %     1:2;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    %     1:2;...
    1:4;...
    1:4;...
    };

NbLayers = [4 6];
NbLayers = NbLayers+1;


% Number of beta images per condition per session (in case each block is
% modeled by a different regressor).
BlockPerCond = 3;

% Options for the SVM
opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions


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
%                            ROIs                           %
% --------------------------------------------------------- %
ROIs_Ori = {...
    'A1', 35, {}, {}; ...
    'TE1.0', 1, {}, {}; ...
    'TE1.1', 2,  {}, {}; ...
    'TE1.2', 3,  {}, {}; ...
    'PT', 4,  {}, {}; ...
    'TE', [1 2 3],  {}, {}; ...
    
    'V1', 11,  {}, {}; ... %5
    'V2', 12,  {}, {}; ...
    'V3d', 13 , {}, {}; ...
    'V3v', 13.5,  {}, {}; ...
    'V4d', 14,  {}, {}; ...
    'V4v', 14.5,  {}, {}; ...
    'V5', 15,  {}, {}; ...
    'V3', [13 13.5],  {}, {}; ...
    'V4', [14 14.5],  {}, {}; ...
    };

opt.verbose = 0;

% --------------------------------------------------------- %
%          Data pre-processing and SVM parameters           %
% --------------------------------------------------------- %
% Normalization and scaling
opt.scaling.img.eucledian = 1;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 1;

% Feature selection (FS)
opt.fs.threshold = 0.75;
opt.fs.type = 'ttest2';

% Recursive feature elminiation (RFE)
opt.rfe.threshold = 0.75;
opt.rfe.nreps = 10;

% SVM C/nu parameters and default arguments
opt.svm.machine = 'C-SVC';
opt.svm.verbose = 1;
if strcmp(opt.svm.machine, 'C-SVC')
    opt.svm.log2c = 2.^(0:9);
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

% Keeps maks and weight output
opt.output.weight = 1;
opt.output.mask = 0;

SaveSufix = '_results_crosslayer';

if opt.fs.do
    SaveSufix = [SaveSufix '_FS'];
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

SaveSufix = [SaveSufix '.mat'];

for SubjInd = [6 7 8 size(SubjectList,1)-1] %1:size(SubjectList,1); %1:size(SubjectList,1)
    
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
        opt.session.nsamples = 3:NbRuns;
    else
        opt.session.nsamples = NbRuns;
    end
    
    AnalysisFolder = fullfile(SubjectFolder, 'Analysis', 'Vertex', 'SVM');
    
    if ~isdir(AnalysisFolder)
        mkdir(AnalysisFolder)
    end
    
    ROIs = ROIs_Ori;
    
    % --------------------------------------------------------- %
    %                     Analysis to perform                   %
    % --------------------------------------------------------- %
%     SVM(1) = struct('name', 'A Stim VS V Stim', 'class', [1 2]);
    
    SVM(1) = struct('name', 'A Stim VS AV Stim', 'class', [1 3]);
    
    SVM(end+1) = struct('name', 'V Stim VS AV Stim', 'class', [2 3]);
    
    SVM(end+1) = struct('name', 'A Att VS V Att', 'class', [5 6]);
    
    SVM(end+1) = struct('name', 'A Stim(A Att VS V Att)', 'class', [7 8]);
    
    SVM(end+1) = struct('name', 'V Stim(A Att VS V Att)', 'class', [9 10]);
    
    SVM(end+1) = struct('name', 'AV Stim(A Att VS V Att)', 'class', [11 12]);
    
    
    %% Get the vertices indices for the ROIs
    
    fprintf(' Reading ROI vertices indices\n')
    
%     delete(fullfile(AnalysisFolder, 'ROIs.mat'))
    
    if exist(fullfile(AnalysisFolder, 'ROIs.mat'), 'file')
        
        load(fullfile(AnalysisFolder, 'ROIs.mat'));
        
    else
        
        SubXP = XP(SubjInd,:);
        
        % LEFT
        % Manually defined A1
        LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
            ['A1_T1_' SubjID '*lcr*']));
        [~, ~, ROI_mapping_A1_L] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
            LogFileList.name), 0, 1); clear LogFileList
        
        ROIs{1,3} = find(ROI_mapping_A1_L==ROIs{1,2}); % Stores the vertices indices
        
        % Cytoarchitechtonicly defined ROIs
        LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'ROI_MNI', ...
            SubXP, [SubXP '-A'], 'SurfaceMeshMapping', ['T1_' SubjID '*lcr*inf_ROI.vtk']));
        [~, ~, ROI_mapping_All_L] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'ROI_MNI', ...
            SubXP, [SubXP '-A'], 'SurfaceMeshMapping', LogFileList.name), 0, 1); clear LogFileList
        
        for iROI = 2:size(ROIs,1)
            temp = [];
            for i=1:length(ROIs{iROI,2})
                temp = [temp find(ROI_mapping_All_L==ROIs{iROI,2}(i))]; % Stores the vertices indices
            end
            ROIs{iROI,3} = temp;
        end
        clear temp
        
        if size(ROI_mapping_All_L,2)~=size(ROI_mapping_A1_L,2)
            error('VTK files for left ROIs have differing number of vertices.')
        end
        
        % Stores the number of vertices in that surface for future sanity check
        NbVertices(1) = size(ROI_mapping_All_L,2);
        
        clear ROI_mapping_All_L ROI_mapping_A1_L
        
        
        % RIGHT
        LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
            ['A1_T1_' SubjID '*rcr*']));
        
        if isempty(LogFileList)
            ROIs{1,4} = [];
        else
            [~, ~, ROI_mapping_A1_R] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'T1_Mapping', ...
                LogFileList.name), 0, 1); clear LogFileList
            
            ROIs{1,4} = find(ROI_mapping_A1_R==ROIs{1,2}+100);
        end
        
        LogFileList = dir(fullfile(SubjectFolder, 'Structural', 'CBS', 'ROI_MNI', ...
            SubXP, [SubXP '-B'], 'SurfaceMeshMapping', ['T1_' SubjID '*rcr*inf_ROI.vtk']));
        [~, ~, ROI_mapping_All_R] = read_vtk(fullfile (SubjectFolder, 'Structural', 'CBS', 'ROI_MNI', ...
            SubXP, [SubXP '-B'], 'SurfaceMeshMapping', LogFileList.name), 0, 1); clear LogFileList
        
        for iROI = 2:size(ROIs,1)
            temp = [];
            for i=1:length(ROIs{iROI,2})
                temp = [temp find(ROI_mapping_All_R==ROIs{iROI,2}(i)+100)];
            end
            ROIs{iROI,4} = temp;
        end
        clear temp
        
        if exist('ROI_mapping_A1_R', 'var')
            if size(ROI_mapping_All_R,2)~=size(ROI_mapping_A1_R,2)
                error('VTK files for right ROIs have differing number of vertices.')
            end
        end
        
        NbVertices(2) = size(ROI_mapping_All_R,2);
        clear ROI_mapping_All_R ROI_mapping_A1_R SubXP
        
        save(fullfile(AnalysisFolder, 'ROIs.mat'), 'NbVertices', 'ROIs');
        
    end
    
    disp(ROIs)
    
    
    %% Creates a cell that lists the names of the beta images as well as their column number in the design matrix
    load(fullfile(SubjectFolder, 'FFX_Block', 'SPM.mat'))
    
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
    
    BetaOfInterest = ~any([BetaNames(:,7)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
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
    CV_Mat_Orig = dataset(...
        {zeros(NbRuns*sum([Class.nbetas]), 1), 'class'}, ...
        {zeros(NbRuns*sum([Class.nbetas]), 1), 'session'});
    
    % For each condition of each class we figure out what is the associated
    % regressors and in which sessions they occur.
    irow = 1; iFile = 1;
    for iClass=1:numel(Class)
        for iCond=1:numel(Class(iClass).cond)
            
            tmp=BetaNames(BetaOfInterest,1:length(char(Class(iClass).cond(iCond))));
            TEMP = BetaOfInterest(strcmp(Class(iClass).cond(iCond), cellstr(tmp)));
            
            for i=1:length(TEMP)
                CV_Mat_Orig.class(irow) = iClass;
                [I,~] = find(TEMP(i)==RegNumbers);
                CV_Mat_Orig.session(irow) = I;
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
    
%     delete(fullfile(AnalysisFolder, 'AllFeatures.mat'))
    
    if exist(fullfile(AnalysisFolder, 'AllFeatures.mat'), 'file')
        
        load(fullfile(AnalysisFolder, 'AllFeatures.mat'));
        
    else
        
        % List which beta image has been mapped onto which vtk file
        [Betas, FilesCell] = CorresBetaVTK(fullfile(SubjectFolder,'Structural', ...
            'CBS','Vertex','fMRI','Beta'));
        
        % Initialize variables
        FeaturesAll = cell(size(ROIs,1),2,numel(NbLayers));
        for iLayer = 1:numel(NbLayers)
            for hs = 1:2
                for iROI=1:size(ROIs,1)
                    FeaturesAll{iROI,hs,iLayer} = nan(size(BetaList,1), numel(ROIs{iROI,hs+2})*NbLayers(iLayer));
                end
            end
        end
        
%         % Set up a waitbar
        Counter = 0;
        t=[];
        h = waitbar(0,'Reading features...','Name', 'Reading features',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
        Steps = numel(NbLayers)*2*size(Betas,1);
        
        % Gets the value for each ROI from the mapping of each beta image
        for iLayer = 1:numel(NbLayers) % the mapping depends on the number of layers
            
            % Format for reading the vertices from the VTK file
            Spec = repmat('%f ', 1, NbLayers(iLayer));
            
            % For the 2 hemispheres
            for hs = 1:2
                
                for iBeta = 1:size(Betas,1)
                    
                    tic;
                    
                    A = fileread(FilesCell{iBeta,hs,iLayer}); % reads file quickly
                    B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
                    
                    C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
                    mapping = cell2mat(C); %clear C
                    
                    if size(mapping,1)~=(NbVertices(hs))
                        error('A VTK file has wrong number of vertices:\n%s', FilesCell{iBeta,hs,iLayer})
                    end
                    
                    % This loop distribute the mapping value to each appropriate ROI
                    for iROI=1:size(ROIs,1)
                        
                        TEMP = mapping(ROIs{iROI,hs+2},:);
                        
                        % This 2 following line are there to remove form further
                        % analysis any vertex that would have at least one 0 or NaN
                        % along its normal.
                        TEMP(TEMP==0)=NaN;
                        TEMP(any(isnan(TEMP),2),:)=nan(sum(any(isnan(TEMP),2)),NbLayers(iLayer));
                        
                        TEMP = TEMP';
                        TEMP = TEMP(:);
                        
                        if all(isnan(TEMP))
                            TEMP = ones(size(TEMP))*Inf;
                        end
                        
                        FeaturesAll{iROI,hs,iLayer}(strcmp(BetaList,Betas(iBeta)),:) = ...
                            repmat(TEMP', sum(strcmp(BetaList,Betas(iBeta))), 1);
                        
                        % clear TEMP
                        
                    end
                    % clear mapping
                    
                    % Update waitbar
                    Counter=Counter+1;
                    if getappdata(h,'canceling')
                        delete(h)
                        return
                    end
                    t(end+1) = toc;
                    waitbar(Counter / Steps, h, sprintf('Subject %s ; %2.2f percent done ; ETA %2.2f min', ...
                        SubjID, Counter/Steps*100, (Steps-Counter)*mean(t)/60));
                    
                    % F = findall(0,'type','figure','tag','TMWWaitbar');
                    % delete(F);
                    
                end
                
            end
            
        end
        
        delete(h)
        
        save(fullfile(AnalysisFolder, 'AllFeatures.mat'), 'FeaturesAll', '-v7.3');
        
    end
    
    for iLayer = 1:numel(NbLayers)
        for hs = 1:2
            for iROI=1:size(ROIs,1)
                
                if ~all(FeaturesAll{iROI,hs,iLayer}(:)==Inf)
                    
                    if ~isempty(ROIs{iROI,hs+2}) && any(all(isnan(FeaturesAll{iROI,hs,iLayer}),2))
                        error('Data for beta image %i is missing.\n', ...
                            str2double(unique(BetaList(all(isnan(FeaturesAll{iROI,hs,iLayer}),2)))))
                    end
                    
                end
                
            end
        end
    end
    
    clear iBeta iROI iLayer hs iLayer t Counter Steps h Spec Betas
    
    delete(fullfile(AnalysisFolder, ['SVM_*' SaveSufix]))
    
    
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
                
                fprintf('   Running on for %i layers\n', NbLayers(iLayer))
                
                FeaturesL = FeaturesAll{iROI,1,iLayer};
                LogFeatL = any(~isnan(FeaturesL));
                FeaturesLayersL = repmat(NbLayers(iLayer):-1:1, 1, size(FeaturesL,2)/NbLayers(iLayer));
                
                FeaturesR = FeaturesAll{iROI,2,iLayer};
                LogFeatR = any(~isnan(FeaturesR));
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
                    
                    if NbSess2Incl < NbRuns
                        fprintf('    Running learning curve with %i sessions\n', NbSess2Incl)
                    else
                        fprintf('    Running analysis with all sessions\n')
                    end
                    
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
                        else
                            NbPerm = opt.permutation.nreps;
                        end
                        
                        %%
                        for iPerm=1:NbPerm
                            
                            CV_Mat = CV_Mat_Orig;
                            
                            %% Permute class within sessions when all sessions are included
                            if iPerm > 1
                                for iRun=1:max(CV_Mat.session)
                                    temp = CV_Mat.class(all([ismember(CV_Mat.class, SVM(iSVM).class), ...
                                        ismember(CV_Mat.session, iRun)], 2));
                                    
                                    CV_Mat.class(all([ismember(CV_Mat.class, SVM(iSVM).class), ...
                                        ismember(CV_Mat.session, iRun)],2)) = temp(randperm(length(temp)));
                                end
                            end
                            temp = [];
                            
                            %% Leave-one-run-out (LORO) cross-validation
                            for iCV=1:size(CV_id, 2)
                                
                                fprintf('     Cross-Validation %i of %i\n', iCV, size(CV_id, 2))
                                
                                TestSess = []; %#ok<NASGU>
                                TrainSess = []; %#ok<NASGU>
                                % Separate training and test sessions
                                [TestSess, TrainSess] = deal(false(size(1:NbRuns)));
                                TestSess(CV_id(iSubSampSess,iCV)) = 1;
                                TrainSess(setdiff(CV_id(iSubSampSess,:), CV_id(iSubSampSess,iCV))) = 1;
                                
                                % Run SVM Left
                                %                                 fprintf('      Left\n')
                                [SVM, weight{1}(iCV,:,:)] = RunSVM(SVM, 1, iLayer, FeaturesL, LogFeatL, FeaturesLayersL, ...
                                    CV_Mat, TrainSess, TestSess, opt, iSVM, iROI, NbSess2Incl, iSubSampSess, iPerm, iCV);
                                
                                % Run SVM Right
                                %                                 fprintf('      Right\n')
                                [SVM, weight{2}(iCV,:,:)] = RunSVM(SVM, 2, iLayer, FeaturesR, LogFeatR, FeaturesLayersR, ...
                                    CV_Mat, TrainSess, TestSess, opt, iSVM, iROI, NbSess2Incl, iSubSampSess, iPerm, iCV);
                                
                                % Run SVM Both
                                %                                 fprintf('      Both hemisphere\n')
                                [SVM, weight{3}(iCV,:,:)] = RunSVM(SVM, 3, iLayer, FeaturesBoth, LogFeatBoth, FeaturesLayersBoth, ...
                                    CV_Mat, TrainSess, TestSess, opt, iSVM, iROI, NbSess2Incl, iSubSampSess, iPerm, iCV);
                                
                                
                                %                                 fprintf('\n')
                                
                            end
                            
                            % Compute weight mean std and sem for each
                            % layer
                            FeaturesLayers{1} = FeaturesLayersL;
                            FeaturesLayers{2} = FeaturesLayersR;
                            FeaturesLayers{3} = FeaturesLayersBoth;
                            
                            LogFeat{1} = LogFeatL;
                            LogFeat{2} = LogFeatR;
                            LogFeat{3} = LogFeatBoth;
                            
                            for iSide = 1:3
                                
                                tmp = squeeze(nanmean(weight{iSide},1));
                                
                                if ~isnan(tmp)
                                    
                                    for i=1:NbLayers(iLayer)
                                        
                                        tmp2 = FeaturesLayers{iSide}(LogFeat{iSide})==i;
                                        
                                        % For the cross layer MVPA
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.layers.mean(i) = nanmean(tmp(tmp2,1));
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.layers.sem(i) = nansem(tmp(tmp2,1));
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.layers.std(i) = nanstd(tmp(tmp2,1));
                                        
                                        % For the whole ROI
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.mean(i) = nanmean(tmp(tmp2,2));
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.sem(i) = nansem(tmp(tmp2,2));
                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl). ...
                                            rand(iSubSampSess).perm(iPerm).weight.std(i) = nanstd(tmp(tmp2,2));
                                        
                                        % clear tmp2
                                    end
                                    
                                end
                                %clear i  tmp
                                
                            end
                            clear iCV TestSess TrainSess results acc weight iSide
                            
                        end % iPerm=1:NbPerm
                        %clear iPerm
                        
                    end % iSubSampSess=1:size(CV_id, 1)
                    %clear iSubSampSess
                    
                end % NbSess2Incl = opt.session.nsamples
                %clear NbSess2Incl
                
            end % iLayer = 1:numel(NbLayers)
            %clear iLayer
            
            %% Calculate prediction accuracies
            for iLayer = 1:numel(NbLayers)
                for Side = 1:3
                    Class_Acc.TotAcc(Side,iLayer) = ...
                        nanmean([SVM(iSVM).ROI(iROI).session(end).rand.perm.CV(:,Side,iLayer).acc]);
                    for iCV=1:size(CV_id, 2)
                        temp(:,:,iCV) = SVM(iSVM).ROI(iROI).session(end).rand.perm.CV(iCV,Side,iLayer).layers.acc;
                    end
                    Class_Acc.TotAccLayers{Side,iLayer} = nanmean(temp,3);
                    temp = [];
                end
                
                if opt.verbose
                    % Display some results
                    fprintf('\n   Accuracy %i layers\n', NbLayers(iLayer));
                    disp(Class_Acc.TotAcc(:,iLayer))
                    disp(Class_Acc.TotAccLayers{1,iLayer})
                    disp(Class_Acc.TotAccLayers{2,iLayer})
                    disp(Class_Acc.TotAccLayers{3,iLayer})
                end
                
            end
            
            % Save data
            Results = SVM(iSVM).ROI(iROI);
            SaveResults(AnalysisFolder, Results, opt, Class_Acc, SVM, iSVM, iROI, SaveSufix)
            
        end % iSubROI=1:numel(SVM(iSVM).ROI)
        clear Mask Features
        
    end % iSVM=1:numel(SVM)
    clear iSVM SVM FeaturesAll
    
end

end

function SaveResults(SaveDir, Results, opt, Class_Acc, SVM, iSVM, iROI, SaveSufix) %#ok<INUSL>

save(fullfile(SaveDir, ['SVM_' SVM(iSVM).name '_ROI_' SVM(iSVM).ROI(iROI).name SaveSufix]), 'Results', 'opt', 'Class_Acc', '-v7.3');

end

function [SVM, weight] = RunSVM(SVM, Side, iLayer, Features, LogFeat, FeaturesLayers, CV_Mat, TrainSess, TestSess, opt, iSVM, iROI, NbSess2Incl, iSubSampSess, iPerm, iCV)

if isempty(Features) || all(Features(:)==Inf)
    
    warning('Empty ROI')
    
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side,iLayer).acc = NaN;
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side, iLayer).layers.acc = NaN;
    
    weight = [];
    
else
    
    [acc, weight, results] = machine_SVC_layers(SVM(iSVM), ...
        Features(:,LogFeat), FeaturesLayers(:,LogFeat), CV_Mat, TrainSess, TestSess, opt);
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side, iLayer).layers.results = {results};
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side, iLayer).layers.acc = acc;
    
    if opt.verbose
        fprintf('\n       Running on all layers.')
    end
    
    results = machine_SVC(SVM(iSVM), Features(:,LogFeat), CV_Mat, TrainSess, TestSess, opt);
    
    weight(1,:,2) = [results.w / sum(results.w)]'; %#ok<NBRAK>
    
    results = rmfield(results, 'w');
    
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side,iLayer).results = {results};
    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,Side,iLayer).acc = mean(results.pred==results.label);
    
end

end
