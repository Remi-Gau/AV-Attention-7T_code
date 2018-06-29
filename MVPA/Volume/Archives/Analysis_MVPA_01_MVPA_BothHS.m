function Analysis_MVPA
clc; clear;

StartDirectory = pwd;

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %         '14';...
    %         '15';...
    '16'
    ];

% Analysing subject 04
% Running SVM:  Auditory areas: V Stim(A Att VS V Att)

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
    %     1:4;...
    1:4;...
    };


% Number of beta images per condition per session (in case each block is
% modeled by a different regressor).
BlockPerCond = 3;

% Options for the SVM
opt.fs.do = 1; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

% Saves an image with...
% ...the voxel kept for the final model
opt.output.mask = 0;
% ...the weight of each of those voxel
opt.output.weight = 1;

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
ROIs = { ...
    'A_lh'; ...
    'A_rh'; ...
    'V_lh'; ...
    'V_rh'; ...
    };


SubROIs = {...
    'V1'; ...
    'V2'; ...
    'V3_4_d'; ...
    'V3_4_v'; ...
    'V5'; ...
    
    'A1'; ...
    'TE'; ...
    'PT'};


HS_Sufix = {...
    '_lh','_L_MNI.nii';
    '_rh','_R_MNI.nii'};

% --------------------------------------------------------- %
%          Data pre-processing and SVM parameters           %
% --------------------------------------------------------- %
% Normalization and scaling
opt.scaling.img.eucledian = 1;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 1;

% Feature selection (FS)
opt.fs.threshold = 0.65;
opt.fs.type = 'ttest2';

% Recursive feature elminiation (RFE)
opt.rfe.threshold = 0.75;
opt.rfe.nreps = 10;

% SVM C/nu parameters and default arguments
opt.svm.machine = 'C-SVC';
if strcmp(opt.svm.machine, 'C-SVC')
    opt.svm.log2c = 0:15;
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


for Restrict = 0:1
    if Restrict==0
        FigSufix='';
    else
        FigSufix='_Restrict';
    end
    
    % Options for the layers
    opt.roi.restrict = Restrict; % Restrict to voxels that are cortical voxels according to CSB tools
    
    SaveSufix = '_results';
    
    if opt.roi.restrict
        SaveSufix = [SaveSufix '_Restrict'];
    end
    
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
    
    for SubjInd = 1%:size(SubjectList,1)
        
        
        % --------------------------------------------------------- %
        %                     Analysis to perform                   %
        % --------------------------------------------------------- %
%         SVM(1) = struct('name', 'Visual areas: A Stim VS V Stim', 'class', [1 2], 'ROI_Nb',  1:5, ...
%             'dir', 'V');
        SVM(1) = struct('name', 'Auditory areas: A Stim VS V Stim', 'class', [1 2], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        
        SVM(end+1) = struct('name', 'Visual areas: A Stim VS AV Stim', 'class', [1 3], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: A Stim VS AV Stim', 'class', [1 3], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        SVM(end+1) = struct('name', 'Visual areas: V Stim VS AV Stim', 'class', [2 3], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: V Stim VS AV Stim', 'class', [2 3], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        SVM(end+1) = struct('name', 'Visual areas: A Att VS V Att', 'class', [5 6], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: A Att VS V Att', 'class', [5 6], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        SVM(end+1) = struct('name', 'Visual areas: A Stim(A Att VS V Att)', 'class', [7 8], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: A Stim(A Att VS V Att)', 'class', [7 8], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        SVM(end+1) = struct('name', 'Visual areas: V Stim(A Att VS V Att)', 'class', [9 10], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: V Stim(A Att VS V Att)', 'class', [9 10], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        SVM(end+1) = struct('name', 'Visual areas: AV Stim(A Att VS V Att)', 'class', [11 12], 'ROI_Nb',  1:5, ...
            'dir', 'V');
        SVM(end+1) = struct('name', 'Auditory areas: AV Stim(A Att VS V Att)', 'class', [11 12], 'ROI_Nb',  6:8, ...
            'dir', 'A');
        
        
        
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
        
        
        %% Creates a cell that lists the full names of the different beta images
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
        
        
        %% Run cross-validation for each model and ROI
        for iSVM=1:numel(SVM)
            
            AnalysisFolder = fullfile(SubjectFolder, 'Analysis', 'ROI', SVM(iSVM).dir);
            
            [~,~,~]=mkdir(AnalysisFolder);
            
            SaveDir = fullfile(AnalysisFolder, 'SVM');
            
            [~,~,~]=mkdir(SaveDir);
            
            %% Creates a dataset that lists for each beta of interest:
            %   - its corresponding class
            %   - the session in which it occurs
            CV_Mat_Orig = dataset(...
                {zeros(NbRuns*sum([Class.nbetas]), 1), 'class'}, ...
                {zeros(NbRuns*sum([Class.nbetas]), 1), 'session'});
            
            FilesList = {};
            
            for hs = 1:2
                
                % For each condition of each class we figure out what is the associated
                % regressors and in which sessions they occur.
                irow = 1;
                
                iFile = 1;
                
                for iClass=1:numel(Class)
                    for iCond=1:numel(Class(iClass).cond)
                        
                        tmp=BetaNames(BetaOfInterest,1:length(char(Class(iClass).cond(iCond))));
                        TEMP = BetaOfInterest(strcmp(Class(iClass).cond(iCond), cellstr(tmp)));
                        
                        for i=1:length(TEMP)
                            CV_Mat_Orig.class(irow) = iClass;
                            [I,~] = find(TEMP(i)==RegNumbers);
                            CV_Mat_Orig.session(irow) = I;
                            irow = irow + 1;
                            
                            FilesList{iFile,hs} = fullfile([AnalysisFolder HS_Sufix{hs,1}], ...
                                sprintf('rbeta_%04d.nii', TEMP(i))); %#ok<*SAGROW>
                            iFile = iFile + 1;
                        end
                        clear TEMP I i tmp
                        
                    end
                end
                
                clear irow iClass iCond
                
                
                %% Gets global mask from GLM
                Mask(hs).global.hdr = spm_vol(fullfile([AnalysisFolder HS_Sufix{hs,1}], 'mask.nii'));
                Mask(hs).global.img = logical(spm_read_vols(Mask(hs).global.hdr));
                
                if opt.roi.restrict
                    tmp = logical(spm_read_vols(spm_vol(fullfile([AnalysisFolder HS_Sufix{hs,1}], 'T1_04_Layers.nii'))));
                    tmp(:,:,:,2) = Mask(hs).global.img;
                    tmp = all(tmp,4);
                    Mask(hs).global.img = tmp;
                    clear tmp
                end
                
                for i=1:length(SVM(iSVM).ROI_Nb)
                    SVM(iSVM).ROI(i).name = SubROIs{SVM(iSVM).ROI_Nb(i),1};
                    Mask(hs).ROI(i).name = SubROIs{SVM(iSVM).ROI_Nb(i),1};
                    Mask(hs).ROI(i).fname = ['rsub_' SubROIs{SVM(iSVM).ROI_Nb(i),1} HS_Sufix{hs,2}];
                    if  SubROIs{SVM(iSVM).ROI_Nb(i),1}(1:2)=='A1' %#ok<STCMP>
                        Mask(hs).ROI(i).fname(end-7:end-4) = [];
                    end
                    Mask(hs).ROI(i).hdr = spm_vol(fullfile([AnalysisFolder HS_Sufix{hs,1}], [Mask(hs).ROI(i).fname]));
                    
                end
                
                hdr = cat(1, Mask(hs).ROI.hdr);
                sts = spm_check_orientations([Mask(hs).global.hdr; hdr]);
                if sts ~= 1
                    error('Images not in same space!');
                end
                
                clear sts hdr i
                
                
                %% Create mask in XYZ format (both world and voxel coordinates)
                [X, Y, Z] = ind2sub(size(Mask(hs).global.img), find(Mask(hs).global.img));
                Mask(hs).global.XYZ = [X'; Y'; Z']; % XYZ format
                clear X Y Z
                Mask(hs).global.size = size(Mask(hs).global.XYZ, 2);
                
                
                %% Combine masks
                for i=1:length(Mask(hs).ROI)
                    tmp = Mask(hs).global.img;
                    tmp(:,:,:,2) = spm_read_vols(Mask(hs).ROI(i).hdr);
                    tmp = all(tmp,4);
                    [X, Y, Z] = ind2sub(size(tmp), find(tmp));
                    
                    Mask(hs).ROI(i).XYZ = [X'; Y'; Z'];
                    Mask(hs).ROI(i).size = size(Mask(hs).ROI(i).XYZ, 2);
                    clear X Y Z tmp
                end
                clear i
                
                for j=1:numel(Mask(hs).ROI)
                    SVM(iSVM).ROI(j).size(hs) = Mask(hs).ROI(j).size;
                end
                clear j
                
                %% Get layers files
                LayerFiles{1,hs} = fullfile([AnalysisFolder HS_Sufix{hs,1}], 'T1_04_Layers.nii'); %#ok<*SAGROW>
                LayerFiles{2,hs} = fullfile([AnalysisFolder HS_Sufix{hs,1}], 'T1_06_Layers.nii');
                LayerFiles{3,hs} = fullfile([AnalysisFolder HS_Sufix{hs,1}], 'T1_10_Layers.nii');
                
                
            end
            
            fprintf('\n Running SVM:  %s\n', SVM(iSVM).name)
            
            %%
            for iSubROI=1:numel(SVM(iSVM).ROI)
                
                %% Read features
                
                % Left HS
                V = spm_vol(char(FilesList{:,1}));
                tmp = spm_get_data(V, Mask(1).ROI(iSubROI).XYZ);
                
                % Right HS
                V = spm_vol(char(FilesList{:,2}));
                tmp2 = spm_get_data(V, Mask(2).ROI(iSubROI).XYZ);
                
                Features = [tmp tmp2];
                clear tmp tmp2 V
                
                %% Get layer labels
                
                % Left HS
                V = spm_vol(char(LayerFiles{:,1}));
                tmp = spm_get_data(V, Mask(1).ROI(iSubROI).XYZ);
                
                % Right HS
                V = spm_vol(char(LayerFiles{:,2}));
                tmp2 = spm_get_data(V, Mask(2).ROI(iSubROI).XYZ);
                
                Layers = [tmp tmp2];
                clear tmp tmp2 V
                
                Layers(:,any(isnan(Features))) = [];
                Features(:,any(isnan(Features))) = [];

                
                % RNG init
                rng('default');
                opt.seed = rng;
                
                Class_Acc = struct('Pred', [], 'Label', [], 'Acc', [], 'TotAcc', []);
                
                
                fprintf('\n  Running ROI:  %s\n', SVM(iSVM).ROI(iSubROI).name)
                fprintf('  Number of voxel before FS/RFE: %i\n', sum(SVM(iSVM).ROI(iSubROI).size))
                
                
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
                            
                            %% Permute class within sessions when all sessions are included
                            if iPerm > 1
                                for iRun=1:max(CV_Mat.session)
                                    temp = CV_Mat.class(all([ismember(CV_Mat.class, SVM(iSVM).class), ...
                                        ismember(CV_Mat.session, iRun)], 2));
                                    
                                    CV_Mat.class(all([ismember(CV_Mat.class, SVM(iSVM).class), ...
                                        ismember(CV_Mat.session, iRun)],2)) = temp(randperm(length(temp)));
                                end
                            end
                            clear temp
                            
                            %% Leave-one-run-out (LORO) cross-validation
                            for iCV=1:size(CV_id, 2)
                                TestSess = [];
                                TrainSess = [];
                                % Separate training and test sessions
                                [TestSess, TrainSess] = deal(false(size(1:NbRuns)));
                                TestSess(CV_id(iSubSampSess,iCV)) = 1;
                                TrainSess(setdiff(CV_id(iSubSampSess,:), CV_id(iSubSampSess,iCV))) = 1;
                                
                                % Run SVM
                                results = machine_SVC(SVM(iSVM), Features, CV_Mat, TrainSess, TestSess, opt);
                                
                                %Saves an image of the weight
                                if opt.output.weight;
                                    
                                    Weights(iCV,:) = [results.w / sum(results.w)]';
                                    
                                    results = rmfield(results, 'w');
                                end
                                
                                SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(iCV) = results;
                                
                                clear results
                                
                            end
                            clear iCV TestSess TrainSess
                            
                            %% Calculate prediction accuracies
                            pred = cat(2, SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV.pred);
                            label = cat(2, SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV.label);
                            
                            Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,iPerm) = mean(mean(pred==label));
                            
                            
                            %% Display some results
                            if iPerm == 1 && NbSess2Incl == NbRuns
                                if ~isempty(SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).fs)
                                    fprintf('  Number of voxel after FS :  %i\n', ...
                                        SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).fs.size);
                                end
                                if  ~isempty(SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).rfe)
                                    fprintf('  Number of voxel after RFE:  %i\n', ...
                                        SVM(iSVM).ROI(iSubROI).session(NbSess2Incl).rand(iSubSampSess).perm(iPerm).CV(1).rfe.size);
                                end
                                fprintf('  Accuracy = %2.3f\n', mean(mean(pred==label)));
                                
                                if opt.output.weight;
                                    
                                    [WeightMean, WeightStd, WeightSEM] =  PlotWeights(SaveDir,SVM(iSVM).name,SVM(iSVM).ROI(iSubROI).name,Weights,Layers,'off',FigSufix);
                                   
                                    SVM(iSVM).ROI(iSubROI).WeightMean=WeightMean;
                                    SVM(iSVM).ROI(iSubROI).WeightStd=WeightStd;
                                    SVM(iSVM).ROI(iSubROI).WeightSEM=WeightSEM;
                                    
                                end
                                
                                if opt.permutation.test
                                    fprintf('  Running permutations\n')
                                end
                            end
                            
                            if opt.permutation.test && iPerm == opt.permutation.nreps && NbSess2Incl == NbRuns
                                fprintf('  Accuracy over permutations = %2.3f +/- %2.3f\n', ...
                                    mean(squeeze(Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,2:iPerm))), ...
                                    std(squeeze(Class_Acc.TotAcc(NbSess2Incl,iSubSampSess,2:iPerm))) )
                            end
                            
                            clear pred label Weights
                            
                        end % iPerm=1:NbPerm
                        
                        clear iPerm
                        
                    end % iSubSampSess=1:size(CV_id, 1)
                    
                    clear iSubSampSess
                    
                end % NbSess2Incl = opt.session.nsamples
                
                clear NbSess2Incl
                
                % Save data into partially readable mat file (see matfile)
                Results = SVM(iSVM).ROI(iSubROI);
                SaveResults(SaveDir, Results, opt, Class_Acc, SVM, iSVM, iSubROI, SaveSufix)
                
            end % iSubROI=1:numel(SVM(iSVM).ROI)
            
            clear Mask Features Layers
            
        end % iSVM=1:numel(SVM)
        
        clear iSVM SVM
        
    end
    
end

end

function SaveResults(SaveDir, Results, opt, Class_Acc, SVM, iSVM, iROI, SaveSufix) %#ok<INUSL>

save(fullfile(SaveDir, ['SVM_' SVM(iSVM).name '_ROI_' SVM(iSVM).ROI(iROI).name SaveSufix]), 'Results', 'opt', 'Class_Acc', '-v7.3');

end


function [WeightMean, WeightStd, WeightSEM] = PlotWeights(SaveDir,SVM,ROI,Weights,Layers,Visible,FigSufix)

Weights = nanmean(Weights);
Weights(isnan(Weights)) = 0;


figure('Name', [SVM ' / ' ROI ' : weights distribution'], ...
    'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible); %#ok<UNRCH>

iSubplot = 1;

FontSize = 10;

for iLayerVol = 1:size(Layers,1)
    
    LayerVol=Layers(iLayerVol,:);
    NbLayers = max(LayerVol(:));
    
    Values2Plot={};
    
    for iLayer=NbLayers:-1:0
        
        tmp = LayerVol==iLayer;
        tmp(2,:) = logical(Weights);
        tmp = all(tmp,1);
        
        WeightMean{iLayerVol,1}(iLayer+1) = mean(Weights(tmp)); %#ok<*NASGU>
        WeightStd{iLayerVol,1}(iLayer+1) = std(Weights(tmp));
        WeightSEM{iLayerVol,1}(iLayer+1) = nansem(Weights(tmp));
        
        Values2Plot{end+1,1} = Weights(tmp); %#ok<*UNRCH>
        
    end
    clear iLayer
    
    subplot(size(Layers,1),size(Weights,4),iSubplot)
    
    grid on
    hold on
    
    distributionPlot(Values2Plot, 'showMM', 1, 'histOpt', 0, 'divFactor', 70, ...
        'globalNorm', 0, 'color', 'b')
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers+1 ,'xticklabel', NbLayers:-1:0, ...
        'ticklength', [0.01 0.01], 'fontsize', FontSize)
    
    t=ylabel('Weight');
    set(t,'fontsize', FontSize);
    
    t=xlabel('Layers');
    set(t,'fontsize', FontSize);
    
    clear Values2Plot
    
    iSubplot = iSubplot+1;
    
end
clear iLayerVol


print(gcf, fullfile(SaveDir, [SVM '_ROI_' ROI FigSufix '.tif']), '-dtiff')


end