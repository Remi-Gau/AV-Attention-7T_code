function Analysis_MVPA_01_MVPA_AVvsAandV
clc; clear;

StartDirectory = pwd;

addpath(fullfile(StartDirectory, 'SubFun'))

SubjectList = [...
    %     '02';...
    '03';...
    %     '04';...
    %'06';...
    '07';...
    '08';...
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
    1:4;...
    %     1:4;...
    %1:2;...
        1:4;...
    1:4;...
    %     1:4;...
    1:4;...
    %     1:4;...
    %     1:4;...
    %1:2;...
    %     1:4;...
    %     1:4;...
    };

FFX = {'0'; '3'; '6'};

NbLayers = 6;

NbWorkers = 4;

% Number of beta images per condition per session (in case each block is
% modeled by a different regressor).
BlockPerCond = 3;

% Options for the SVM
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


% Class(4) = struct('name', 'A+V Stim', 'cond', cell(1), 'nbetas', 2*BlockPerCond);
% Class(end).cond = {'A+V Stim - Auditory Attention' 'A+V Stim - Visual Attention'};

% 
% Class(5) = struct('name', 'A Att', 'cond', cell(1), 'nbetas', 3*BlockPerCond);
% Class(end).cond = {'A Stim - Auditory Attention' 'V Stim - Auditory Attention' 'AV Stim - Auditory Attention'};
% 
% Class(6) = struct('name', 'V Att', 'cond', cell(1), 'nbetas', 3*BlockPerCond);
% Class(end).cond = {'A Stim - Visual Attention' 'V Stim - Visual Attention' 'AV Stim - Visual Attention'};
% 
% 
% 
% Class(7) = struct('name', 'A Att A Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'A Stim - Auditory Attention'};
% 
% Class(8) = struct('name', 'V Att A Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'A Stim - Visual Attention'};
% 
% 
% Class(9) = struct('name', 'A Att V Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'V Stim - Auditory Attention'};
% 
% Class(10) = struct('name', 'V Att V Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'V Stim - Visual Attention'};
% 
% 
% Class(11) = struct('name', 'A Att AV Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'AV Stim - Auditory Attention'};
% 
% Class(12) = struct('name', 'V Att AV Stim', 'cond', cell(1), 'nbetas', 1*BlockPerCond);
% Class(end).cond = {'AV Stim - Visual Attention'};

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
    opt.svm.log2c = 1;
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
            
            for Norm = 6
                
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
                
                
                for iFFX = 3 %1:length(FFX)
                    
                    SaveSufix = '_results_vol_FixedC';
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
                    
                    SaveSufix = [SaveSufix '_FWHM_' FFX{iFFX} '_Layers_' num2str(NbLayers) '.mat'];
                    
                    for SubjInd = 1:size(SubjectList,1)
                        
                        % --------------------------------------------------------- %
                        %                            ROIs                           %
                        % --------------------------------------------------------- %
                        
%                         Mask.ROI(1) = struct('name', 'TE1.0', 'fname', 'rrwTE_1.0_MNI.nii');
%                         Mask.ROI(end+1) = struct('name', 'TE1.1', 'fname', 'rrwTE_1.1_MNI.nii');
%                         Mask.ROI(end+1) = struct('name', 'TE1.2', 'fname', 'rrwTE_1.2_MNI.nii');
                        Mask.ROI(1) = struct('name', 'TE', 'fname', 'rrwTE_MNI.nii');
%                         Mask.ROI(end+1) = struct('name', 'STG', 'fname',  'rrwHG_STG_AAL.nii');
                        
                        Mask.ROI(end+1) = struct('name', 'V1', 'fname', 'rrwProbRet_V1.nii');
                        Mask.ROI(end+1) = struct('name', 'V2', 'fname', 'rrwProbRet_V2.nii');
                        Mask.ROI(end+1) = struct('name', 'V3', 'fname', 'rrwProbRet_V3.nii');
%                         Mask.ROI(end+1) = struct('name', 'V4', 'fname', 'rrwProbRet_V4.nii');
%                         Mask.ROI(end+1) = struct('name', 'V5', 'fname', 'rrwProbRet_V_hMT.nii');
                        
                        Mask.ROI(end+1) = struct('name', 'STGpost', 'fname', 'rrwSTG_Post_AAL.nii');
                        Mask.ROI(end+1) = struct('name', 'V1-2-3', 'fname', 'rrwProbRet_V1-2-3.nii');
                        
                        
                        % --------------------------------------------------------- %
                        %                     Analysis to perform                   %
                        % --------------------------------------------------------- %
                        SVM(1) = struct('name', 'AV VS A+V', 'class', [3 104], 'ROI',  1:length(Mask.ROI), 'CombineClass', 1:2);
                        
                        
                        % --------------------------------------------------------- %
                        %                        Subject data                       %
                        % --------------------------------------------------------- %
                        SubjID = SubjectList(SubjInd,:);
                        
                        fprintf('\n\nAnalysing subject %s\n', SubjID)
                        
                        SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
                        
                        RoiFolder = fullfile(SubjectFolder, 'Transfer', 'ROI');
                        
                        % AnalysisFolder = fullfile(SubjectFolder, 'FFX', ['AVT_HD_' FFX{iFFX}]);
                        AnalysisFolder = fullfile(SubjectFolder, 'Transfer');
                        
                        SaveDir = fullfile(AnalysisFolder, 'SVM');
                        if ~isdir(SaveDir)
                            mkdir(SaveDir)
                        end
                        
                        NbRuns = length(FoldersNames{SubjInd});
                        
                        % If we want to have a learning curve
                        if opt.session.curve
                            % #sessions over which to run the learning curve
                            opt.session.nsamples = 3:NbRuns;
                        else
                            opt.session.nsamples = NbRuns;
                        end
                        
                        
                        %% Gets global mask from GLM and ROI masks for the data
                        Mask.global.hdr = spm_vol(fullfile(AnalysisFolder, 'rmask.nii'));
                        Mask.global.img = logical(spm_read_vols(Mask.global.hdr));
                        
                        for i=1:length(Mask.ROI)
                            Mask.ROI(i).hdr = spm_vol(fullfile(RoiFolder, Mask.ROI(i).fname));
                        end
                        
                        hdr = cat(1, Mask.ROI.hdr);
                        sts = spm_check_orientations([Mask.global.hdr; hdr]);
                        if sts ~= 1
                            error('Images not in same space!');
                        end
                        
                        clear sts hdr i
                        
                        
                        %% Create mask in XYZ format (both world and voxel coordinates)
                        [X, Y, Z] = ind2sub(size(Mask.global.img), find(Mask.global.img));
                        Mask.global.XYZ = [X'; Y'; Z']; % XYZ format
                        clear X Y Z
                        Mask.global.size = size(Mask.global.XYZ, 2);
                        Mask.global.XYZmm = Mask.global.hdr.mat(1:3,:) ...
                            * [Mask.global.XYZ; ones(1, Mask.global.size)]; % voxel to world transformation
                        
                        
                        %% Combine masks
                        xY.def = 'mask';
                        for i=1:length(Mask.ROI)
                            xY.spec = fullfile(RoiFolder, Mask.ROI(i).fname);
                            [xY, Mask.ROI(i).XYZmm, j] = spm_ROI(xY, Mask.global.XYZmm);
                            Mask.ROI(i).XYZ = Mask.global.XYZ(:,j);
                            Mask.ROI(i).size = size(Mask.ROI(i).XYZ, 2);
                        end
                        
                        clear xY j i
                        
                        
                        %% Gets Layer labels
                        fprintf(' Reading layer labels\n')
                        
                        LayerLabelsFile = fullfile(AnalysisFolder, ['T1_0' num2str(NbLayers) '_Layers.nii']);
                        
                        LayerLabelsHdr = spm_vol(LayerLabelsFile);
                        
                        sts = spm_check_orientations([Mask.global.hdr; LayerLabelsHdr]);
                        if sts ~= 1
                            error('Images not in same space!');
                        end
                        clear sts
                        
                        
                        for i=1:length(Mask.ROI)
                            LayerLabels{i} = spm_get_data(LayerLabelsHdr, Mask.ROI(i).XYZ);
                        end
                        
                        
                        %% Creates a cell that lists the full names of the different beta images
                        load(fullfile(AnalysisFolder, 'SPM.mat'))
                        
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
                        
                        FilesList = {};
                        
                        % For each condition of each class we figure out what is the associated
                        % regressors and in which sessions they occur.
                        
                        if iFFX==1
                            PreFix=[];
                        elseif iFFX==2
                            PreFix='S3';
                        elseif iFFX==3
                            PreFix='S6';
                        end
                        
                        irow = 1;
                        for iClass=1:numel(Class)
                            for iCond=1:numel(Class(iClass).cond)
                                
                                tmp=BetaNames(BetaOfInterest,1:length(char(Class(iClass).cond(iCond))));
                                TEMP = BetaOfInterest(strcmp(Class(iClass).cond(iCond), cellstr(tmp)));
                                
                                for i=1:length(TEMP)
                                    CV_Mat_Orig(irow,1) = iClass;
                                    [I,~] = find(TEMP(i)==RegNumbers);
                                    CV_Mat_Orig(irow,2) = I;
                                    irow = irow + 1;
                                    
                                    FilesList{end+1,1} = fullfile(AnalysisFolder, ...
                                        sprintf('%srbeta_%04d.nii', PreFix, TEMP(i))); %#ok<*SAGROW>
                                end
                                clear TEMP I i tmp
                                
                            end
                        end
                        clear irow iClass iCond
                        
                        % Mask each image by each ROI and create a features set (images x voxel)
                        Features = cell(1, length(Mask.ROI));
                        
                        fprintf(' Reading features\n')
                        V = spm_vol(char(FilesList));
                        
                        parfor i=1:length(Mask.ROI)
                            tmp = spm_get_data(V, Mask.ROI(i).XYZ);
                            %                             tmp(:,any(isnan(tmp)))=[];
                            Features{i} = tmp;
                        end
                        
                        % Keep one hearder in case we need it later
                        HDR = struct('dim', V(1).dim, 'dt', V(1).dt, 'mat', V(1).mat, 'pinfo', V(1).pinfo);
                        
                        clear V irow RegNumbers RegNames iClass iCond
                        
                        
                        %% Run cross-validation for each model and ROI
                        for i=1:numel(SVM)
                            for j=1:numel(SVM(i).ROI)
                                SVM(i).ROI_XYZ{j,1} = Mask.ROI(SVM(i).ROI(j)).XYZ;
                            end
                            SVM(i).ROI = struct('name', {Mask.ROI(SVM(i).ROI).name}, 'size', {Mask.ROI(SVM(i).ROI).size});
                        end
                        
                        clear i j
                        
                        for iSVM=1:numel(SVM)
                            
                            for iROI=1:numel(SVM(iSVM).ROI)
                                
                                % RNG init
                                rng('default');
                                opt.seed = rng;
                                
                                Class_Acc = struct('Pred', [], 'Label', [], 'Acc', [], 'TotAcc', []);
                                
                                % ROI index
                                ROI = ismember({Mask.ROI.name}, SVM(iSVM).ROI(iROI).name);
                                
                                
                                fprintf('\nAnalysing subject %s\n', SubjID)
                                fprintf(' Running SVM:  %s\n', SVM(iSVM).name)
                                fprintf('  Running ROI:  %s\n', SVM(iSVM).ROI(iROI).name)
                                fprintf('  Number of voxel before FS/RFE: %i\n', SVM(iSVM).ROI(iROI).size)
                                
                                
                                %% Subsample sessions for the learning curve (otherwise take
                                % all of them)
                                for NbSess2Incl = opt.session.nsamples
                                    
                                    if NbSess2Incl < NbRuns
                                        fprintf('  Running learning curve with %i sessions\n', NbSess2Incl)
                                    else
                                        fprintf('  Running analysis with all sessions\n')
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
                                            
                                            Features_ROI = Features{ROI};
                                            
                                            
                                            % Combine features if needed
                                            if isfield(SVM(iSVM), 'CombineClass')
                                                
                                                Sessions = unique(CV_Mat(:,2));
                                                Class2Create = SVM(iSVM).class>100;
                                                Class2Fuse =  SVM(iSVM).CombineClass;
                                                
                                                for iSess = 1:numel(Sessions)
                                                    Row2Sum1 = find(all([ismember(CV_Mat(:,1),Class2Fuse(1)) CV_Mat(:,2)==Sessions(iSess)],2));
                                                    Row2Sum2 = find(all([ismember(CV_Mat(:,1),Class2Fuse(2)) CV_Mat(:,2)==Sessions(iSess)],2));
                                                    for i=1:numel(Row2Sum1)
                                                        Features_ROI(end+1,:) = sum(Features_ROI([Row2Sum1(i) Row2Sum2(i)],:));
                                                        CV_Mat(end+1,:) = [SVM(iSVM).class(Class2Create) Sessions(iSess)];
                                                    end
                                                    clear Row2Sum1 Row2Sum2
                                                end

                                                clear tmp Sessions Class2Create Class2Fuse
                                            end
                                            
                                            
                                            
                                            %% Permute class within sessions when all sessions are included
                                            if iPerm > 1
                                                for iRun=1:max(CV_Mat(:,2))
                                                    temp = CV_Mat((all([ismember(CV_Mat(:,1), SVM(iSVM).class), ...
                                                        ismember(CV_Mat(:,2), iRun)], 2)),1);
                                                    
                                                    CV_Mat(:,1) = CV_Mat(temp(randperm(length(temp))),1);
                                                end
                                            end
                                            clear temp
                                            
                                            LayerLabels_ROI = LayerLabels{ROI};

                                            LayerLabels_ROI(:,any(isnan(Features_ROI),1))=[];
                                            Features_ROI(:,any(isnan(Features_ROI),1))=[];
                                            
                                            Features_ROI = Features_ROI(:,logical(LayerLabels_ROI));
                                            LayerLabels_ROI = LayerLabels_ROI(:,logical(LayerLabels_ROI));
                                            
                                            %% Leave-one-run-out (LORO) cross-validation
                                            tic;
                                            for iCV=1:size(CV_id, 2)
                                                TestSess = []; %#ok<NASGU>
                                                TrainSess = []; %#ok<NASGU>
                                                % Separate training and test sessions
                                                [TestSess, TrainSess] = deal(false(size(1:NbRuns)));
                                                TestSess(CV_id(iSubSampSess,iCV)) = 1;
                                                TrainSess(setdiff(CV_id(iSubSampSess,:), CV_id(iSubSampSess,iCV))) = 1;
                                                
                                                % Run SVM
                                                TEMP{iCV,1} = machine_SVC(SVM(iSVM), Features_ROI, CV_Mat, TrainSess, TestSess, opt);
                                                
                                                % Run SVM
                                                for iLayer= 1:NbLayers
                                                    Features_ROI_Layer = Features_ROI(:,LayerLabels_ROI==iLayer);
                                                    TEMP2{iCV,iLayer} = machine_SVC(SVM(iSVM), Features_ROI_Layer, CV_Mat, TrainSess, TestSess, opt);
                                                end
                                            end
                                            clear iCV TestSess TrainSess
                                            
                                            t = toc;
                                            SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).ExecTime = t;
                                            
                                            for iCV=1:size(CV_id, 2)
                                                SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,1) = TEMP{iCV,1};
                                                for iLayer = 1:NbLayers
                                                    SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV,1+iLayer) = TEMP2{iCV,iLayer};
                                                end
                                            end
                                            clear TEMP
                                            
                                            % Calculate prediction accuracies
                                            for iLayer = 1:(NbLayers+1)
                                                pred = cat(2, SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(:,iLayer).pred);
                                                label = cat(2, SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(:,iLayer).label);
                                                Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,iPerm,iLayer) = mean(mean(pred==label));
                                            end
                                            
                                            %% Display some results
                                            if iPerm == 1 && NbSess2Incl == NbRuns
                                                if ~isempty(SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).fs)
                                                    fprintf('  Number of voxel after FS :  %i\n', ...
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).fs.size);
                                                end
                                                if  ~isempty(SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).rfe)
                                                    fprintf('  Number of voxel after RFE:  %i\n', ...
                                                        SVM(iSVM).ROI(iROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).rfe.size);
                                                end
                                                seq = repmat(' %0.3f',[1 NbLayers+1]);
                                                fprintf(['  Accuracy\n' seq '\n'], squeeze(Class_Acc.TotAcc(end,1,1,:))');
                                                if opt.permutation.test
                                                    fprintf('  Running permutations\n')
                                                end
                                            end
                                            
                                            if opt.permutation.test && iPerm == opt.permutation.nreps && NbSess2Incl == NbRuns
                                                fprintf('  Accuracy over permutations = %2.3f +/- %2.3f\n', ...
                                                    mean(squeeze(Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,2:iPerm))), ...
                                                    std(squeeze(Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,2:iPerm))) )
                                            end
                                            
                                            pred =[]; label =[];
                                            
                                        end % iPerm=1:NbPerm
                                        
                                        clear iPerm
                                        
                                    end % iSubSampSess=1:size(CV_id, 1)
                                    
                                    clear iSubSampSess
                                    
                                end % NbSess2Incl = opt.session.nsamples
                                
                                clear NbSess2Incl
                                
                                % Save data into partially readable mat file (see matfile)
                                Results = SVM(iSVM).ROI(iROI);
                                SVM(iSVM).ROI(iROI).session=[]; % Remove the last results to save memory.
                                SaveResults(SaveDir, Results, opt, Class_Acc, SVM, iSVM, iROI, SaveSufix)
                                
                            end % iSubROI=1:numel(SVM(iSVM).ROI)
                            
                        end % iSVM=1:numel(SVM)
                        
                        clear iSVM SVM Mask Features
                        
                    end
                    
                end
                
            end
        end
    end
end


if str2double(MatlabVer(1:4))>2013
    delete(gcp);
else
    matlabpool close
end

end


function SaveResults(SaveDir, Results, opt, Class_Acc, SVM, iSVM, iROI, SaveSufix) %#ok<INUSL>

save(fullfile(SaveDir, ['SVM_' SVM(iSVM).name '_ROI_' SVM(iSVM).ROI(iROI).name SaveSufix]), 'Results', 'opt', 'Class_Acc', '-v7.3');

end
