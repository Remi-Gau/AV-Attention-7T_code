clc; clear; close all;

PlotSubjects = 0;
Switch = 1;

NbLayers = 6;
NbLayersMVPA = 6;


% StartDirectory='/home/rxg243/Dropbox/PhD/Experiments/AV_Integration_7T/Archives';
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




%% Plot activations and deactivations

close all
clear DATA


DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;

figure('position', FigDim, 'name', 'activations_and_deactivations', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);


DATA.MVPA = 0;
DATA.YLabelInset = 1;
DATA.InsetLim = [.5 .25;-.5 -.15];
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
subplot(4,4,2)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,2)
iCond = 2;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));


PlotProfileAndBetas(DATA)

ax = subplot(4,4,6);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute;
PlotInsetFinal(DATA)





% PT
subplot(4,4,10)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,10)
iCond = 2;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,14);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)


DATA.InsetLim = [.5 .1;-2 -.5];
if DATA.WithSubj
    DATA.MIN = -3;
    DATA.MAX = 0.1;
else
    DATA.MIN = -1.6;
    DATA.MAX = 0.1;
end

% V1
subplot(4,4,4)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,4)
iCond = 1;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,8);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)



% V2-3
subplot(4,4,12)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,12)
iCond = 1;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,16);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)




% Plot Activations

DATA.YLabel = 'B Param. est. [a u]';

DATA.InsetLim = [2.7 0.78;-.1 -.02];

if DATA.WithSubj
    DATA.MIN = -0.3;
    DATA.MAX = 4.5;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end

DATA.OneSideTTest = {'right' 'both' 'both'};


% TE
subplot(4,4,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,1)
iCond = 1;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));


PlotProfileAndBetas(DATA)

ax = subplot(4,4,5);
axis('off')
DATA.ax = ax.Position;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'B Param. est. [a u]';

% PT
subplot(4,4,9)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,9)
iCond = 1;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,13);
axis('off')
DATA.ax = ax.Position;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = '';

DATA = rmfield(DATA, 'YLabel');

if DATA.WithSubj
    DATA.MIN = -0.3;
    DATA.MAX = 4.2;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end

% V1
subplot(4,4,3)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,3)
iCond = 2;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,7);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)

% V2-3
subplot(4,4,11)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,11)
iCond = 2;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,15);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)



% print(gcf, fullfile(FigureFolder,['Fig2_activations_and_deactivations_' MedianSufix '.tif']), '-dtiff')

print(gcf, fullfile(FigureFolder,['Fig2_activations_and_deactivations_' MedianSufix '.tif']), '-dtiff')


%% Plot Cross Modal Influence
close all
clear DATA

DATA.WithSubj = 0;%PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;


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
DATA.InsetLim = [1.2 0.3; -1.2 -.3];
if DATA.WithSubj
    DATA.MIN = -1.3;
    DATA.MAX = 1.6;
else
    DATA.MIN = -0.25;
    DATA.MAX = 0.7;
end


% TE
DATA.YLabel = 'B Param. est. [a u]';
subplot(4,4,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,1)
iCond = 3;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,5);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 1;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'B Param. est. [a u]';

% PT
subplot(4,4,9)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,9)
iCond = 3;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,13);
axis('off')
DATA.ax = ax.Position;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'B Param. est. [a u]';

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
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,7);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'B Param. est. [a u]';

% V2-3
subplot(4,4,11)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(4,4,11)
iCond = 4;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,15);
axis('off')
DATA.ax = ax.Position;
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'B Param. est. [a u]';


% Plot MVPA
DATA.MVPA = 1;
DATA.InsetLim = [0.4 0.08; -.4 -.08];
if DATA.WithSubj
    DATA.MIN = 0.25;
    DATA.MAX = 0.95;
else
    DATA.MIN = 0.45;
    DATA.MAX = 0.8;
end

DATA.YLabel = 'Decoding accuracy';

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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';

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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';


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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';

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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)


print(gcf, fullfile(FigureFolder,['Fig3_CrossModalInfluence' MedianSufix '.tif']), '-dtiff')



%% Plot attention effects for type of stimulus
close all

Stim = {'All'};

BOLD_Cdtion = 5:8;
MVPA_Cdtion = 3:6;

for iStim = 1
    
    clear DATA
    
    DATA.WithSubj = PlotSubjects;
    
    DATA.Scatter = Scatter;
    DATA.WithPerm = 0;
    DATA.PlotInset = 0;
    DATA.YLabel = 'B Param. est. [a u]';
    
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
    
    % TE
    subplot(4,4,1)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,1)
    iCond = BOLD_Cdtion(iStim);
    iROI = 1;
    DATA.Name = '';
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,5);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 1;
    DATA.YLabel = 'S Param. est. [a u]';
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    DATA.YLabel = 'B Param. est. [a u]';
    
    
    % PT
    subplot(4,4,9)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,9)
    iCond = BOLD_Cdtion(iStim);
    iROI = 2;
    DATA.Name = '';
    DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,13);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabel = 'S Param. est. [a u]';
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    DATA.YLabel = 'B Param. est. [a u]';
    
    
    % V1
    subplot(4,4,3)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,3)
    iCond = BOLD_Cdtion(iStim);
    iROI = 3;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,7);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.YLabel = 'S Param. est. [a u]';
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    DATA.YLabel = 'B Param. est. [a u]';
    
    
    % V2-3
    subplot(4,4,11)
    PlotRectangle(NbLayers,Fontsize,Switch)
    subplot(4,4,11)
    iCond = BOLD_Cdtion(iStim);
    iROI = 4;
    DATA.Name = '';
    DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
    DATA.Color =  [0 0 0];
    DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,1:2,iCond));
    DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
    PlotProfileAndBetas(DATA)
    
    ax = subplot(4,4,15);
    axis('off')
    DATA.ax = ax.Position;
    DATA.YLabelInset = 0;
    DATA.YLabel = 'S Param. est. [a u]';
    DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
    DATA.YLabel = 'B Param. est. [a u]';
    
    
    % Plot MVPA
    DATA = rmfield(DATA, 'OneSideTTest');
    
    DATA.MVPA = 1;
    DATA.YLabel = 'Decoding accuracy';
    
    DATA.InsetLim = [0.38 .07; -.38 -.07];
    if DATA.WithSubj
        DATA.MIN = .15;
        DATA.MAX = .85;
    else
        DATA.MIN = 0.4;
        DATA.MAX = 0.7;
    end
    
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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';
    
    
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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';
    
    
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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';
    
    
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
DATA.YLabel = 'S Param. est. [a u]';
DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
DATA.YLabel = 'Decoding accuracy';
    
    print(gcf, fullfile(FigureFolder,['Fig4_Attention_Stim_' Stim{iStim} '_' MedianSufix '.tif']), '-dtiff')

    
end

