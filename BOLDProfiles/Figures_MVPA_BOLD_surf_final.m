clc; clear; close all;

PlotSubjects = 0;
Switch = 1;

NbLayers = 6;
NbLayersMVPA = 6;


StartDirectory=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartDirectory, 'SubFun')))
Get_dependencies('/home/rxg243/Dropbox/')
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

FigDim = [100 100 1800 1000];
Visible = 'on';
Transparent = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);


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
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
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
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
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
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
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
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
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
        [~,P] = ttest(tmp);
        All_P(iROI,:,iSVM) = P;
        
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
    Legends1 = {'', 'Constant', '','','', '', '', '', '', 'Linear'};
    Legends2 = {...
        'ROI', ...
        'mean', '(','STD',')', 't value','p value', 'effect size', '', ...
        'mean', '(','STD',')', 't value','p value', 'effect size'};
elseif size(DesMat,2)==3
    Legends1 = {'', 'Constant', '','','', '', '', '', '', 'Linear', '','','', '', '', '', '', 'Quad'};
    Legends2 = {...
        'ROI', ...
        'mean', '(','STD',')', 't value','p value', 'effect size', '', ...
        'mean', '(','STD',')', 't value','p value', 'effect size', '', ...
        'mean', '(','STD',')', 't value','p value', 'effect size'};
end

%% Plot Subjects Legend

figure(1)
hold on
COLOR_Subject = ColorSubject();
for iSubj=1:size(COLOR_Subject,1)
    plot([0 1], [iSubj iSubj], 'color', COLOR_Subject(iSubj,:), 'linewidth', 2)
end
legend(SubjectList)

print(gcf, fullfile(FigureFolder,'SubjectLegend.tif'), '-dtiff')

close all


%% Plot deactivations
SavedTxt = fullfile(FigureFolder,['Deactivations' MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');

DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

close all

figure('position', FigDim, 'name', 'Deactivations', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

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
DATA.YLabelInset = 1;
DATA.InsetLim = [.5 .4;-.5 -.4];
if DATA.WithSubj
    DATA.MIN = -0.6;
    DATA.MAX = 1;
else
    DATA.MIN = -0.3;
    DATA.MAX = 0.15;
end

% DATA.OneSideTTest = {'left' 'left' 'both'};
DATA.OneSideTTest = {'left' 'both' 'both'};


% A1
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,1)
iCond = 2;
iROI = 1;
DATA.Name = char({'V vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));


PlotProfileAndBetas(DATA)

ax = subplot(4,2,3);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute;
PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,5)
iCond = 2;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

DATA.InsetLim = [.5 .4;-2 -.4];
if DATA.WithSubj
    DATA.MIN = -3;
    DATA.MAX = 0.1;
else
    DATA.MIN = -1.6;
    DATA.MAX = 0.1;
end

% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,2)
iCond = 1;
iROI = 3;
DATA.Name = char({'A vs. Baseline';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,4);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(4,2,6)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,6)
iCond = 1;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,8);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

print(gcf, fullfile(FigureFolder,['Fig2_Deactivations' MedianSufix '.tif']), '-dtiff')

fclose (fid);



%% Plot Activations
clear DATA
SavedTxt = fullfile(FigureFolder,['Activations' MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');

DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

close all

figure('position', FigDim, 'name', 'Activations', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

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
DATA.YLabelInset = 1;

DATA.InsetLim = [3.5 0.9;-.1 -.1];

if DATA.WithSubj
    DATA.MIN = -0.3;
    DATA.MAX = 4.5;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end


% DATA.OneSideTTest = {'right' 'left' 'both'};
DATA.OneSideTTest = {'right' 'both' 'both'};


% TE
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,1)
iCond = 1;
iROI = 1;
DATA.Name = char({'A vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));


PlotProfileAndBetas(DATA)

ax = subplot(4,2,3);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,5)
iCond = 1;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

if DATA.WithSubj
    DATA.MIN = -0.3;
    DATA.MAX = 4.2;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end

% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,2)
iCond = 2;
iROI = 3;
DATA.Name = char({'V vs. Baseline';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,4);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(4,2,6)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,2,6)
iCond = 2;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,8);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


print(gcf, fullfile(FigureFolder,['Fig1_Activations' MedianSufix '.tif']), '-dtiff')

fclose (fid);


%% Plot Cross Modal Influence
close all
clear DATA
SavedTxt = fullfile(FigureFolder,['CrossModal' MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');

DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

figure('position', FigDim, 'name', 'Cross Modal Influence', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);


% Plot BOLD
DATA.MVPA = 0;
DATA.InsetLim = [1.1 0.3; -1 -.2];
if DATA.WithSubj
    DATA.MIN = -1.3;
    DATA.MAX = 1.6;
else
    DATA.MIN = -0.25;
    DATA.MAX = 0.7;
end

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
subplot(4,4,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,1)
iCond = 3;
iROI = 1;
DATA.Name = char({'AV vs. A';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,5);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 1;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)


% PT
subplot(4,4,9)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,9)
iCond = 3;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,13);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)

if DATA.WithSubj
    DATA.MIN = -1.2;
    DATA.MAX = 1;
else
    DATA.MIN = -0.35;
    DATA.MAX = 0.1;
end

% V1
subplot(4,4,3)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,3)
iCond = 4;
iROI = 3;
DATA.Name = char({'AV vs. V';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,7);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)


% V2-3
subplot(4,4,11)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,11)
iCond = 4;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,15);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)



% Plot MVPA
DATA.MVPA = 1;
DATA.InsetLim = [0.4 0.07; -.2 -.07];
if DATA.WithSubj
    DATA.MIN = 0.25;
    DATA.MAX = 0.95;
else
    DATA.MIN = 0.45;
    DATA.MAX = 0.8;
end

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
subplot(4,4,2)
PlotRectangle(NbLayersMVPA,Fontsize,Switch)
subplot(4,4,2)
iSVM = 1;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));


PlotProfileAndBetas(DATA)

ax = subplot(4,4,6);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)

% PT
subplot(4,4,10)
PlotRectangle(NbLayersMVPA,Fontsize,Switch)
subplot(4,4,10)
iSVM = 1;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,14);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)


% V1
subplot(4,4,4)
PlotRectangle(NbLayersMVPA,Fontsize,Switch)
subplot(4,4,4)
iSVM = 2;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,8);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)


% V2-3
subplot(4,4,12)
PlotRectangle(NbLayersMVPA,Fontsize,Switch)
subplot(4,4,12)
iSVM = 2;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,16);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

Print2Table(fid, ROIs, iROI, DATA)


print(gcf, fullfile(FigureFolder,['Fig3_CrossModalInfluence' MedianSufix '.tif']), '-dtiff')

fclose (fid);


%% Plot attention effects for type of stimulus
close all

Stim = {'All', 'A', 'V', 'AV'};

BOLD_Cdtion = 5:8;
MVPA_Cdtion = 3:6;

for iStim = 2:3
    
    clear DATA
    SavedTxt = fullfile(FigureFolder,['Attention_Stim_' Stim{iStim} MedianSufix '.csv']);
    fid = fopen (SavedTxt, 'w');
    
    DATA.WithSubj = PlotSubjects;
    
    DATA.Scatter = Scatter;
    DATA.WithPerm = 0;
    DATA.PlotInset = 0;
    DATA.YLabel = 'Param. est. [a u]';
    
    figure('position', FigDim, 'name', 'Attention', 'Color', [1 1 1], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    % Plot BOLD
    DATA.MVPA = 0;
    DATA.InsetLim = [0.8 0.2; -.8 -.2];
    if DATA.WithSubj
        DATA.MIN = -1.5;
        DATA.MAX = 2;
    else
        DATA.MIN = -0.3;
        DATA.MAX = 0.5;
    end
    
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
    subplot(4,4,1)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,1)
    iCond = BOLD_Cdtion(iStim);
    iROI = 1;
    DATA.Name = char({['Stim ' Stim{iStim} ' - A vs. V att'];'';'A1'});
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,5);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 1;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % PT
    subplot(4,4,9)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,9)
    iCond = BOLD_Cdtion(iStim);
    iROI = 2;
    DATA.Name = 'PT';
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,13);
    axis('off')
    DATA.ax = ax.Position;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % V1
    subplot(4,4,3)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,3)
    iCond = BOLD_Cdtion(iStim);
    iROI = 3;
    DATA.Name = char({['Stim ' Stim{iStim} ' - V vs. A att'];'';'V1'});
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,7);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % V2-3
    subplot(4,4,11)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,11)
    iCond = BOLD_Cdtion(iStim);
    iROI = 4;
    DATA.Name = 'V2-3';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,15);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % Plot MVPA
    DATA = rmfield(DATA, 'OneSideTTest');
    
    DATA.MVPA = 1;
    DATA.YLabel = 'Accuracy';
    
    DATA.InsetLim = [0.28 .06; -.28 -.06];
    if DATA.WithSubj
        DATA.MIN = 0;
        DATA.MAX = 1
    else
        DATA.MIN = 0.4;
        DATA.MAX = 0.7;
    end
    
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
    subplot(4,4,2)
    PlotRectangle(NbLayersMVPA,Fontsize,Switch)
    subplot(4,4,2)
    iSVM = MVPA_Cdtion(iStim);
    iROI = 1;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Color =  [0 0 0];
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,6);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA)
    
    
    % PT
    subplot(4,4,10)
    PlotRectangle(NbLayersMVPA,Fontsize,Switch)
    subplot(4,4,10)
    iSVM = MVPA_Cdtion(iStim);
    iROI = 2;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Color =  [0 0 0];
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,14);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA)
    
    
    % V1
    subplot(4,4,4)
    PlotRectangle(NbLayersMVPA,Fontsize,Switch)
    subplot(4,4,4)
    iSVM = MVPA_Cdtion(iStim);
    iROI = 3;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Color =  [0 0 0];
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,8);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA)
    
    
    % V2-3
    subplot(4,4,12)
    PlotRectangle(NbLayersMVPA,Fontsize,Switch)
    subplot(4,4,12)
    iSVM = MVPA_Cdtion(iStim);
    iROI = 4;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,:)));
    DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
    DATA.Color =  [0 0 0];
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,16);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA)
    
    print(gcf, fullfile(FigureFolder,['Fig4_Attention_Stim_' Stim{iStim} '_' MedianSufix '.tif']), '-dtiff')
    
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
    
    DATA.WithSubj = PlotSubjects;
    
    DATA.Scatter = Scatter;
    DATA.WithPerm = 0;
    DATA.PlotInset = 0;
    DATA.YLabel = 'Param. est. [a u]';
    
    figure('position', FigDim, 'name', 'Attention', 'Color', [1 1 1], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    % Plot BOLD
    DATA.MVPA = 0;
    DATA.InsetLim = [1.4 0.3; -.6 -.2];
    if DATA.WithSubj
        DATA.MIN = -2.5;
        DATA.MAX = 2.5;
    else
        DATA.MIN = -.1;
        DATA.MAX = .6;
    end
    
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
    subplot(2,2,1)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(2,2,1)
    iCond = BOLD_Cdtion(iStim);
    iROI = 5;
    DATA.Name = char('V1 act');
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [ 0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    
    PlotProfileAndBetas(DATA)
    
    %     ax = subplot(4,2,3);
    %     axis('off')
%     DATA.ax = ax.Position;
    DATA.YLabelInset = 1;
    DATA.ToPermute = ToPermute;
    % PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    % V2-3 act
    subplot(2,2,2)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(2,2,2)
    iCond = BOLD_Cdtion(iStim);
    iROI = 7;
    DATA.Name = char('V2-3 act');
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [ 0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    %     ax = subplot(4,2,4);
    %     axis('off')
%     DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute;
    % PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    
    % Plot BOLD
    DATA.MVPA = 0;
    DATA.InsetLim = [0.5 0.2; -1.2 -.3];
    if DATA.WithSubj
        DATA.MIN = -2.5;
        DATA.MAX = 2.5;
    else
        DATA.MIN = -.8;
        DATA.MAX = 0.1;
    end
    
    
    % V1 deact
    subplot(2,2,3)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(2,2,3)
    iCond = BOLD_Cdtion(iStim);
    iROI = 6;
    DATA.Name = 'V1 deact';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [ 0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    %     ax = subplot(4,2,7);
    %     axis('off')
%     DATA.ax = ax.Position;
    DATA.ToPermute = ToPermute;
    % PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    % V2-3 deact
    subplot(2,2,4)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(2,2,4)
    iCond = BOLD_Cdtion(iStim);
    iROI = 8;
    DATA.Name = 'V2-3 deact';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [ 0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    %     ax = subplot(4,2,8);
    %     axis('off')
%     DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.ToPermute = ToPermute;
    % PlotInsetFinal(DATA)
    
    Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)
    
    
    print(gcf, fullfile(FigureFolder,['Fig4_Attention_Stim_' Stim{iStim} '_ActDeact' MedianSufix '.tif']), '-dtiff')
    
    fclose (fid);
    
    
end