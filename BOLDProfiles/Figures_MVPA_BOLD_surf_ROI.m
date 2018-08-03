clc; clear; close all;

PlotSubjects = 1;
Switch = 1;

NbLayers = 6;
NbLayersMVPA = 6;


StartDirectory=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartDirectory, 'SubFun')))

% SourceFolder = fullfile('/media/rxg243/BackUp2/AV_Integration_7T_2/Results/Profiles/Surfaces');
SourceFolder = fullfile(StartDirectory, 'Figures', 'ProfilesSurface', strcat(num2str(NbLayers), '_layers'));

FigureFolder = fullfile(StartDirectory, 'Figures', strcat(num2str(NbLayers+2), '_layers'));


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
load(fullfile(SourceFolder, 'Data_BOLD_WholeROI.mat'), 'AllSubjects_Data')


% A against baseline & V against baseline
Target=1;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).MainEffects.DATA(:,iCond);
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
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).BiVSUni.DATA(:,iCond);
    end
    Target = Target+1;
end


% Main effect of Att
Target=5;
for iCond = 2
    for iROI=1:length(AllSubjects_Data)    
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).Differential.MainEffect.DATA(:,iCond);
    end
    Target = Target+1;
end

% Effect of Att for each stim
Target=6;
for iCond = 9:11
    for iROI=1:length(AllSubjects_Data)
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).Differential.DATA(:,iCond);
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

            MVPA_SubjectsData(SubjInd,iSVM,iROI) = mean([CV(:).acc]); %#ok<*SAGROW>
  
            
            clear Acc Pred Label B X
        end
        
    end
    
end


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

Legends2 = {'ROI', 'mean', '','STD','', 'p value', 'effect size'};

%% Plot deactivations
close all

SavedTxt = fullfile(FigureFolder, 'Deactivations_ROI.csv');
fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD \n');
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

figure('position', FigDim, 'name', 'Deactivations', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.MIN = -1.75;
DATA.MAX = 0.25;

DATA.OneSideTTest = 'left';
DATA.ToPermute=ToPermute;

% TE
subplot(2,2,1)
iCond = 2;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'V vs. Baseline';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% PT
subplot(2,2,3)
iCond = 2;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
subplot(2,2,2)
iCond = 1;
iROI = 3;
DATA.Legend{2} = char({'A vs. Baseline';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
subplot(2,2,4)
iCond = 1;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

print(gcf, fullfile(FigureFolder, 'Deactivations_WholeROI.tif'), '-dtiff')

fclose (fid);


%% Plot Activations
clear DATA

SavedTxt = fullfile(FigureFolder,'Activations_ROI.csv');
fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD \n');
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

figure('position', FigDim, 'name', 'Activations', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.MIN = -0.1;
DATA.MAX = 3.2;

DATA.OneSideTTest = {'right'};
DATA.ToPermute=ToPermute;

% TE
subplot(2,2,1)
iCond = 1;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'V vs. Baseline';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% STG
subplot(2,2,3)
iCond = 1;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
subplot(2,2,2)
iCond = 2;
iROI = 3;
DATA.Legend{2} =  char({'A vs. Baseline';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
subplot(2,2,4)
iCond = 2;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


print(gcf, fullfile(FigureFolder, 'Fig2_Activations_WholeROI.tif'), '-dtiff')


fclose (fid);

%% Plot Cross Modal Influence
clear DATA

SavedTxt = fullfile(FigureFolder, 'CrossModal_ROI.csv');
fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD \n');
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

figure('position', FigDim, 'name', 'Cross Modal Influence', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);


DATA.OneSideTTest = {'both'};
DATA.ToPermute=ToPermute;

% Plot BOLD
DATA.MVPA = 0;
DATA.MIN = -1;
DATA.MAX = 0.8;

% TE
subplot(2,4,1)
iCond = 5;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'AV VS A';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% PT
subplot(2,4,5)
iCond = 5;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)


% V1
subplot(2,4,3)
iCond = 6;
iROI = 3;
DATA.Legend{2} = char({'AV VS A';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,7)
iCond = 6;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)



% Plot MVPA
DATA.MIN = 0.4;
DATA.MAX = 0.85;
DATA.MVPA = 1;
DATA.Legend{1} = 'Accuracy';

fprintf (fid, '\n\n');
fprintf (fid, 'MVPA profile\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

% TE
subplot(2,4,2)
iSVM = 1;
iROI = 1;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% PT
subplot(2,4,6)
iSVM = 1;
iROI = 2;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V1
subplot(2,4,4)
iSVM = 2;
iROI = 3;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,8)
iSVM = 2;
iROI = 4;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)


print(gcf, fullfile(FigureFolder,'Fig3_CrossModalInfluence_WholeROI.tif'),  '-dtiff')


fclose (fid);


%% Plot attention effects
clear DATA

SavedTxt = fullfile(FigureFolder,'Attention_ROI.csv');
fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD \n');
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

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
DATA.MIN = -0.8;
DATA.MAX = 0.8;

DATA.OneSideTTest = {'both'};
DATA.ToPermute=ToPermute;


% TE
subplot(2,4,1)
iCond = 7;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'A vs. V att';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI)*-1;
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% STG
subplot(2,4,5)
iCond = 7;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI)*-1;
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V1
subplot(2,4,3)
iCond = 7;
iROI = 3;
DATA.Legend{2} = char({'A vs. V att';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(2,4,7)
iCond = 7;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% Plot MVPA
DATA.MIN = 0.4;
DATA.MAX = 0.85;
DATA.MVPA = 1;

DATA.OneSideTTest = {'both'};
DATA.ToPermute=ToPermute;


fprintf (fid, '\n\n');
fprintf (fid, 'MVPA profile\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

% TE
subplot(2,4,2)
iSVM = 3;
iROI = 1;
DATA.Legend{1} = 'Accuracy';
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% STG
subplot(2,4,6)
iSVM = 3;
iROI = 2;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V1
subplot(2,4,4)
iSVM = 3;
iROI = 3;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,8)
iSVM = 3;
iROI = 4;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
% PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

print(gcf, fullfile(FigureFolder,'Fig4_Attention_WholeROI.tif'), '-dtiff')

fclose (fid);