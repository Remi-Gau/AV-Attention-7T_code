clc; clear; close all;

StartDirectory = fullfile(pwd, '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))


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


ROIs= {...
    'A1_surf'
    'PT_BT'
    'V1_surf'
    'V2-3_surf'
    };

Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');


% Holms Thresholds
HolmsThresholds = .05*ones(1,length(ROIs))./(length(ROIs)+1-(1:length(ROIs)));

FigDim = [100 100 1800 1000];
Visible = 'on';
Transparent = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);
NbLayers = 6;

DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));


%% Get data for BOLD
% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 8 Interaction AV-A under A att vs. V under V att
% 9 Interaction AV-V under A att vs. V under V att
% 10 Interaction between A attention vs. V attention and A vs. V
% load(fullfile(FigureFolder, strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')
load(fullfile(FigureFolder, strcat('Data_BOLD_WholeROI.mat')), 'AllSubjects_Data')

% A against baseline & A against baseline
Target=1;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).MainEffects.DATA(:,iCond);
    end
    Target = Target+1;
end

% Interactions
Target=8;
for iCond = 1:3
    for iROI=1:length(AllSubjects_Data)
        BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).Interaction.DATA(:,iCond);
    end
    Target = Target+1;
end


% 3 A under A att vs. V under V att
% 4 V under A att vs. V under V att
% 5 AV-A irrespective of attention
% 6 AV-V irrespective of attention
% 7 A attention vs. V attention (pooled over all stimulation conditions)

% A_A vs A_V   &   V_A vs V_V
Target=3;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)

            BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).Differential.DATA(:,iCond);
    end
    Target = Target+1;
end

% AV-A and AV-V
Target=5;
for iCond = 7:8
    for iROI=1:length(AllSubjects_Data)
            BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).BiVSUni.DATA(:,iCond);
    end
    Target = Target+1;
end

% Main effect of Att
Target=7;
for iCond = 2
    for iROI=1:length(AllSubjects_Data)
            BOLD_SubjectsData(:,Target,iROI) = AllSubjects_Data(iROI).Differential.MainEffect.DATA(:,iCond);
    end
end


%% Get Data for MVPA
cd(FigureFolder)

opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;
opt.scaling.idpdt = 1;

Include = repmat(logical(ones(size(SubjectList,1),1)),[1,numel(ROIs)]);



for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    %%
    for iSVM=1:numel(Analysis)
        
        for iROI=1:numel(ROIs)
            
            Save_vol = [...
                'SVM_' Analysis(iSVM).name...
                '_ROI_' ROIs{iROI}];
            
            SaveSufix = '_results_vol_FixedC';
            
            if opt.scaling.idpdt
                SaveSufix = [SaveSufix '_Idpdt']; %#ok<*AGROW>
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
            
            SaveSufix = [SaveSufix '_FWHM_0Slab_Layers_' num2str(NbLayers) '.mat'];
            
            Save_vol = [Save_vol SaveSufix];
            
            try 
                load(fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID],  'Transfer', 'SVM', ...
                    Save_vol), 'Class_Acc', 'Results');
                MVPA_SubjectsData(SubjInd,iSVM,iROI) = squeeze(Class_Acc.TotAcc(end,:,:,1)); %#ok<*SAGROW>

            catch
                MVPA_SubjectsData(SubjInd,iSVM,iROI) = NaN;

            end
            
            clear Acc Pred Label B X
        end
        
    end
    
end


% BOLD
% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 3 A under A att vs. A under V att
% 4 V under A att vs. V under V att
% 5 AV-A irrespective of attention
% 6 AV-V irrespective of attention
% 7 A attention vs V attention (pooled over all stimulation conditions)
% 8 Interaction between AV vs A for A att vs V att
% 9 Interaction between AV vs V for A att vs V att
% 10 Interaction between A attention vs V attention and A vs V

% MVPA
% Analysis(1) = struct('name', 'A Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');

Legends2 = {'ROI', 'mean', '(','SEM',')', 'p value', 'effect size'};

%% Plot deactivations
close all

SavedTxt = fullfile(StartDirectory, 'Paper','Deactivations_ROI.csv');
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

% TE
subplot(2,2,1)
iCond = 2;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'V vs. Baseline';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% PT
subplot(2,2,3)
iCond = 2;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
subplot(2,2,2)
iCond = 1;
iROI = 3;
DATA.Legend{2} = char({'A vs. Baseline';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
subplot(2,2,4)
iCond = 1;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
     'Deactivations_WholeROI.tif'), '-dtiff')
 
 fclose (fid);

 
%% Plot Activations
clear DATA

SavedTxt = fullfile(StartDirectory, 'Paper','Activations_ROI.csv');
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

% TE
subplot(2,2,1)
iCond = 1;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'V vs. Baseline';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% STG
subplot(2,2,3)
iCond = 1;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V1
subplot(2,2,2)
iCond = 2;
iROI = 3;
DATA.Legend{2} =  char({'A vs. Baseline';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

% V2-3
subplot(2,2,4)
iCond = 2;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
     'Fig2_Activations_WholeROI.tif'), '-dtiff')

 
  fclose (fid);

%% Plot Cross Modal Influence
clear DATA

SavedTxt = fullfile(StartDirectory, 'Paper','CrossModal_ROI.csv');
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
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% PT
subplot(2,4,5)
iCond = 5;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)


% V1
subplot(2,4,3)
iCond = 6;
iROI = 3;
DATA.Legend{2} = char({'AV VS A';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,7)
iCond = 6;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

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
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% PT
subplot(2,4,6)
iSVM = 1;
iROI = 2;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V1
subplot(2,4,4)
iSVM = 2;
iROI = 3;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,8)
iSVM = 2;
iROI = 4;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
     'Fig3_CrossModalInfluence_WholeROI.tif'),  '-dtiff')


  fclose (fid);


%% Plot attention effects
clear DATA

SavedTxt = fullfile(StartDirectory, 'Paper','Attention_ROI.csv');
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
DATA.OneSideTTest = {'right'};


% TE
subplot(2,4,1)
iCond = 7;
iROI = 1;
DATA.Legend{1} = 'Param. est. [a u]';
DATA.Legend{2} = char({'A vs. V att';'';'A1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI)*-1;
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% STG
subplot(2,4,5)
iCond = 7;
iROI = 2;
DATA.Legend{2} = 'PT';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI)*-1;
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V1
subplot(2,4,3)
iCond = 7;
iROI = 3;
DATA.Legend{2} = char({'A vs. V att';'';'V1'});
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(2,4,7)
iCond = 7;
iROI = 4;
DATA.Legend{2} = 'V2-3';
DATA.Data = BOLD_SubjectsData(:,iCond,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% Plot MVPA
DATA.MIN = 0.4;
DATA.MAX = 0.85;
DATA.MVPA = 1;
DATA = rmfield(DATA, 'OneSideTTest');

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
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% STG
subplot(2,4,6)
iSVM = 3;
iROI = 2;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V1
subplot(2,4,4)
iSVM = 3;
iROI = 3;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

% V2-3
subplot(2,4,8)
iSVM = 3;
iROI = 4;
DATA.Legend{2} = '';
DATA.Data = MVPA_SubjectsData(:,iSVM,iROI);
PlotROIForFig(DATA)

Print2TableROI(fid, ROIs, iROI, DATA)

cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
     'Fig4_Attention_WholeROI.tif'), '-dtiff')

   fclose (fid);