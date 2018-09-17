clc; clear; close all;

NbLayers = 6;
NbLayersMVPA = 6;

StartDirectory=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartDirectory, 'AV-Attention-7T_code', 'SubFun')))
% Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

SourceFolder = fullfile(StartDirectory, 'Figures', 'ProfilesSurface', strcat(num2str(NbLayers), '_layers'));

FigureFolder = fullfile(StartDirectory, 'Figures', strcat(num2str(NbLayers+2), '_layers'));
mkdir(FigureFolder)

Median = 1;

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

ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3';...
    'V1_act';...
    'V1_deact';...
    'V23_act';...
    'V23_deact';...
    };

for iSubj=1:size(SubjectList,1)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


if Median
    MedianSufix = '_Median'; %#ok<*UNRCH>
else
    MedianSufix = '';
end

%% Get data for BOLD
% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
load(fullfile(SourceFolder, strcat('Data_Surf', MedianSufix ,'_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A against baseline & V against baseline
Target=1;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA(:,iCond,Include))';
%         [~,P] = ttest(tmp);
%         All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).MainEffects.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% 3 AV-A irrespective of attention
% 4 AV-V irrespective of attention
% 5 A attention vs. V attention (pooled over all stimulation conditions)
% 6 A attention vs. V attention for A stim
% 7 A attention vs. V attention for V stim
% 8 A attention vs. V attention for AV stim

% AV-A and AV-V
Target=3;
for iCond = 7:8
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.DATA(:,iCond,Include))';
%         [~,P] = ttest(tmp);
%         All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).BiVSUni.DATA(:,iCond,iSubj);
            NormBOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).NormBiVSUni.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% Main effect of Att
Target=5;
for iCond = 2
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA(:,iCond,Include))';
%         [~,P] = ttest(tmp);
%         All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = ...
                AllSubjects_Data(iROI).Differential.MainEffect.DATA(1:NbLayers,iCond,iSubj);
            NormBOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = ...
                AllSubjects_Data(iROI).NormDifferential.MainEffect.DATA(1:NbLayers,iCond,iSubj);
        end
    end
end


% Effect of Att for each stim
Target=6;
for iCond = 9:11
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA(:,iCond,Include))';
%         [~,P] = ttest(tmp);
%         All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Differential.DATA(1:NbLayers,iCond,iSubj);
            NormBOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).NormDifferential.DATA(1:NbLayers,iCond,iSubj);
        end
    end
    Target = Target+1;
end

Include =[];
for iROI=1:numel(ROIs)
    if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
        temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
        Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
    end
end
clear AllSubjects_Data temp iROI


%% Get Data for MVPA
cd(SourceFolder)

Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');

opt.acroslayer.do = 0;
opt.leave2out.do = 0;

opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;
opt.scaling.idpdt = 1;

Include = repmat(logical(ones(size(SubjectList,1),1)),[1,numel(ROIs)]);

DesMat = (1:NbLayersMVPA)-mean(1:NbLayersMVPA);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayersMVPA,1) DesMat'];
DesMat = spm_orth(DesMat);


FFX = {'0'};

SaveSufix = '_results_surf_FixedC';
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

SaveSufix = [SaveSufix '_FWHM_' FFX{1} '_Layers_' num2str(NbLayersMVPA+2) '.mat'];


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    %%
    for iSVM=1:numel(Analysis)
        
        for iROI=1:numel(ROIs)
            
            Save_vol = [...
                'SVM_' Analysis(iSVM).name...
                '_ROI_' ROIs{iROI} SaveSufix];
            
            load(fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID],  'Transfer', 'SVM', Save_vol));
            
            CV = Results.session(end).rand.perm.CV;
            
            for iLayer= 2:NbLayersMVPA+1
                Label = [];
                Pred = [];
                for iCV=1:numel(CV)
                    Label(:,end+1) = CV(iCV).layers.results{1}{iLayer}.label;
                    Pred(:,end+1) = CV(iCV).layers.results{1}{iLayer}.pred(:,iLayer);
                end
                MVPA_SubjectsData(SubjInd,iSVM,iROI,iLayer-1) = mean(mean(Pred==Label,2));
                Acc(iLayer-1,:) = mean(Pred==Label,2)';
            end
            
            X=repmat(DesMat,size(Acc,2),1);
            Acc=flipud(Acc(:)-.5);
            [B,~,~] = glmfit(X, Acc, 'normal', 'constant', 'off');
            SubjectsBetas(SubjInd,iSVM,iROI,1:size(X,2)) = B;
            
            clear Acc Pred Label B X
        end
        
    end
    
end

% Compile p values and betas
for iSVM = 1:numel(Analysis)
    for iROI=1:numel(ROIs)
        tmp = squeeze(SubjectsBetas(logical(Include(:,iROI)),iSVM,iROI,:));
%         [~,P] = ttest(tmp);
%         All_P(iROI,:,iSVM) = P;
        
        MVPA_SubjectsBetas(:,iROI,1:size(tmp,2),iSVM) = tmp;
    end
end
clear P SubjectsBetas


% BOLD
% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 3 AV-A irrespective of attention
% 4 AV-V irrespective of attention
% 5 A attention vs. V attention (pooled over all stimulation conditions)
% 6 A attention vs. V attention for A stim
% 7 A attention vs. V attention for V stim
% 8 A attention vs. V attention for AV stim

% MVPA
% Analysis(1) = struct('name', 'A Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');

if size(DesMat,2)==2
    Legends1 = {'', 'Constant', '','','', '', '', '', '',  '',  '',  '', 'Linear'};
    Legends2 = {...
        'ROI', ...
        'mean', '(','STD',')', 't value','p value', 'effect size', '',  '',  '',  '',...
        'mean', '(','STD',')', 't value','p value', 'effect size'};
elseif size(DesMat,2)==3
    Legends1 = {'', 'Constant', '','','', '', '', '', '',  '',  '',  '', 'Linear', '','','', '', '', '', '',  '',  '',  '', 'Quad'};
    Legends2 = {...
        'ROI', ...
        'mean', '(','STD',')', 't value','p value', 'effect size', '',  '',  '',  '',...
        'mean', '(','STD',')', 't value','p value', 'effect size', '',  '',  '',  '',...
        'mean', '(','STD',')', 't value','p value', 'effect size'};
end



%% Plot deactivations
SavedTxt = fullfile(FigureFolder,['Deactivations' MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');

DATA.WithPerm = 0;

fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

DATA.MVPA = 0;
% DATA.OneSideTTest = {'left' 'left' 'both'};
DATA.OneSideTTest = {'left' 'both' 'both'};

% A1
iCond = 2;
iROI = 1;
DATA.Name = char({'V vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% PT
iCond = 2;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
iCond = 1;
iROI = 3;
DATA.Name = char({'A vs. Baseline';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
iCond = 1;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

fclose (fid);



%% Plot Activations
clear DATA
SavedTxt = fullfile(FigureFolder,['Activations' MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');
DATA.WithPerm = 0;

fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

DATA.MVPA = 0;

% DATA.OneSideTTest = {'right' 'left' 'both'};
DATA.OneSideTTest = {'right' 'both' 'both'};

% A1
iCond = 1;
iROI = 1;
DATA.Name = char({'A vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% PT
iCond = 1;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
iCond = 2;
iROI = 3;
DATA.Name = char({'V vs. Baseline';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
iCond = 2;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

fclose (fid);


%% Plot Cross Modal Influence
clear DATA

SavedTxt = fullfile(FigureFolder,['CrossModal' MedianSufix '.csv']);
fid = fopen (SavedTxt, 'w');

DATA.WithPerm = 0;

DATA.MVPA = 0;
fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

% TE
iCond = 3;
iROI = 1;
DATA.Name = char({'AV vs. A';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% PT
iCond = 3;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% V1
iCond = 4;
iROI = 3;
DATA.Name = char({'AV vs. V';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% V2-3
iCond = 4;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)


% Plot MVPA
DATA.MVPA = 1;
DATA.YLabel = 'Decoding accuracy';
fprintf (fid, '\n\n');
fprintf (fid, 'MVPA profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

% TE
iSVM = 1;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% PT
iSVM = 1;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% V1
iSVM = 2;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

% V2-3
iSVM = 2;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
DATA.ToPermute = ToPermute;
Print2Table(fid, ROIs, iROI, DATA)

fclose (fid);


%% Plot attention effects for type of stimulus
close all

Stim = {'All', 'A', 'V', 'AV'};

BOLD_Cdtion = 5:8;
MVPA_Cdtion = 3:6;

for iStim = 1
    
    clear DATA
    
    SavedTxt = fullfile(FigureFolder,['Attention_Stim_' Stim{iStim} MedianSufix '.csv']);
    fid = fopen (SavedTxt, 'w');
    
    DATA.WithPerm = 0;
    
    
    % DATA.OneSideTTest = {'right' 'left' 'both'};
    DATA.OneSideTTest = {'both' 'both' 'both'};

    fprintf (fid, 'BOLD profile\n');
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    % TE
    iCond = BOLD_Cdtion(iStim);
    iROI = 1;
    DATA.Name = char({['Stim ' Stim{iStim} ' - A vs. V att'];'';'A1'});
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % PT
    iCond = BOLD_Cdtion(iStim);
    iROI = 2;
    DATA.Name = 'PT';
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % V1
    iCond = BOLD_Cdtion(iStim);
    iROI = 3;
    DATA.Name = char({['Stim ' Stim{iStim} ' - V vs. A att'];'';'V1'});
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % V2-3
    iCond = BOLD_Cdtion(iStim);
    iROI = 4;
    DATA.Name = 'V2-3';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % Plot MVPA
    DATA = rmfield(DATA, 'OneSideTTest');
    
    DATA.MVPA = 1;
    DATA.YLabel = 'Accuracy';
    
    fprintf (fid, '\n\n');
    fprintf (fid, 'MVPA profile\n');
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    % TE
    iSVM = MVPA_Cdtion(iStim);
    iROI = 1;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA)
    
    
    % PT
    iSVM = MVPA_Cdtion(iStim);
    iROI = 2;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA)
    
    % V1
    iSVM = MVPA_Cdtion(iStim);
    iROI = 3;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA)
    
    % V2-3
    iSVM = MVPA_Cdtion(iStim);
    iROI = 4;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA)

    fclose (fid);
    
end


%% Plot attention effects for type of stimulus act/deact
close all

Stim = {'All', 'A', 'V', 'AV'};

BOLD_Cdtion = 5:8;
MVPA_Cdtion = 3:6;

for iStim = 2:3
    
    clear DATA
    SavedTxt = fullfile(FigureFolder,['Attention_Stim_' Stim{iStim} '_ActDeact' MedianSufix '.csv']);
    fid = fopen (SavedTxt, 'w');
    
    DATA.WithPerm = 0;

    % Plot BOLD
    DATA.MVPA = 0;
    DATA.OneSideTTest = {'both' 'both' 'both'};
    
    fprintf (fid, 'BOLD profile\n');
    for i=1:length(Legends1)
        fprintf (fid, '%s,', Legends1{i});
    end
    fprintf (fid, '\n');
    for i=1:length(Legends2)
        fprintf (fid, '%s,', Legends2{i});
    end
    fprintf (fid, '\n');
    
    % V1 Act
    iCond = BOLD_Cdtion(iStim);
    iROI = 5;
    DATA.Name = char('V1 act');
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % V2-3 act
    iCond = BOLD_Cdtion(iStim);
    iROI = 7;
    DATA.Name = char('V2-3 act');
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % Plot BOLD
    DATA.MVPA = 0;

    % V1 deact
    iCond = BOLD_Cdtion(iStim);
    iROI = 6;
    DATA.Name = 'V1 deact';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % V2-3 deact
    iCond = BOLD_Cdtion(iStim);
    iROI = 8;
    DATA.Name = 'V2-3 deact';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.ToPermute = ToPermute;
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

    fclose (fid);
    
end