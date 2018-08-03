clc; clear; close all;

NbLayers = 6;

StartDirectory=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartDirectory, 'SubFun')))

% StartDirectory = '/home/rxg243/Dropbox/PhD/Experiments/AV_Integration_7T';

SourceFolder = fullfile(StartDirectory, 'Figures', 'Profiles', 'Volumes', strcat(num2str(NbLayers), '_layers'));
% SourceFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));

FigureFolder = fullfile('/data/AV_Integration_2', 'Figures', strcat(num2str(NbLayers), '_layers'));
mkdir(FigureFolder)


ANTs = 0;
Median = 1;

% SubjectList = [...
%     '02';...
%     '03';...
%     '04';...
%     '07';...
%     '08';...
%     '09';...
%     '11';...
%     '12';...
%     '13';...
%     '15';...
%     '16'
%     ];

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


ROIs= {...
    'A1_surf'
    'PT_surf_thres'
    'V1_surf_thres'
    'V2-3_surf_thres'
    };

% ROIs= {...
%     'A1_surf'
%     'PT_BT'
%     'V1_surf'
%     'V2-3_surf'
%     };


% Holms Thresholds
HolmsThresholds = .05*ones(1,length(ROIs))./(length(ROIs)+1-(1:length(ROIs)));

FigDim = [100 100 1800 1000];
Visible = 'on';
Transparent = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);


if Median
    MedianSufix = '_Median'; %#ok<*UNRCH>
else
    MedianSufix = '';
end

if ANTs
    ANTsSufix = '_ANTs';
else
    ANTsSufix = '';
end

%% Get data for BOLD

% 1 A attention vs. V attention for A stim
% 2 A attention vs. V attention for V stim
% 3 A attention vs. V attention for AV stim
load(fullfile(SourceFolder, strcat('Data', MedianSufix ,'_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', ANTsSufix, '.mat')), 'AllSubjects_Data')

% A_A vs A_V   &   V_A vs V_V
Target=1;
for iCond = 9:11
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Differential.DATA(1:NbLayers,iCond,iSubj);
        end
    end
    Target = Target+1;
end

Include =[];
% AllSubjects_Data(2).name='STGpost';
for iROI=1:numel(ROIs)
    if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
        temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
        Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
    end
end
Include([4 11],:) = [];
clear AllSubjects_Data temp iROI

% Holms Thresholds
for iCond = 1:size(All_P,3)
    for iP=1:size(All_P,2)
        tmp = All_P(:,iP,iCond);
        [~,I] = ismember(tmp,sort(tmp,'ascend'));
        BOLD_Thresholds(:,iP,iCond) = HolmsThresholds(I); %#ok<*SAGROW>
    end
end
clear All_P I iP iCond


% BOLD
% 1 A attention vs. V attention for A stim
% 2 A attention vs. V attention for V stim
% 3 A attention vs. V attention for AV stim


Legends1 = {'', 'Constant', '','','', '', '', '', 'Linear'};
Legends2 = {'ROI', 'mean', '(','SEM',')', 'p value', 'effect size', '', 'mean', '(','SEM',')', 'p value', 'effect size'};


%% Plot Attention effect for A stim
DATA.Scatter = Scatter;
DATA.WithSubj = 1;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

SavedTxt = fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_A' ANTsSufix MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

close all

figure('position', FigDim, 'name', 'MainEffect_Att_Stim_A', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.MIN = -1.2;
DATA.MAX = 2.5;
DATA.InsetLim = [1.6 0.3];
DATA.YLabelInset = 1;

DATA.OneSideTTest = {'right' 'left'};


% A1
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,1)
iCond = 1;
iROI = 1;
DATA.Name = char({'A Stim - A vs. V att';'';'A1'});
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Color =  [1 0 0];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(BOLD_SubjectsBetas(:,1,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,2,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,3,:,1)),...
    squeeze(BOLD_SubjectsBetas(:,4,:,1)));

PlotProfileAndBetas(DATA)

ax = subplot(4,2,3); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,5)
iCond = 1;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,2)
iCond = 1;
iROI = 3;
DATA.Name = char({'A Stim - A vs. V att';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,4); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% V2-3
subplot(4,2,6)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,6)
iCond = 1;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0.5 0.5 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,8); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_A' ANTsSufix MedianSufix '.tif']), '-dtiff')
 
 
fclose (fid);


%% Plot Attention effect for V stim
DATA.Scatter = Scatter;
DATA.WithSubj = 1;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

close all

SavedTxt = fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_V' ANTsSufix MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

figure('position', FigDim, 'name', 'MainEffect_Att_Stim_V', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.MIN = -.9;
DATA.MAX = .9;
DATA.InsetLim = [1 0.2];
DATA.YLabelInset = 1;

DATA.OneSideTTest = {'right' 'left'};


% A1
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,1)
iCond = 2;
iROI = 1;
DATA.Name = char({'V Stim - A vs. V att';'';'A1'});
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Color =  [1 0 0];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(BOLD_SubjectsBetas(:,1,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,2,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,3,:,1)),...
    squeeze(BOLD_SubjectsBetas(:,4,:,1)));

PlotProfileAndBetas(DATA)

ax = subplot(4,2,3); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,5)
iCond = 2;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,2)
iCond = 2;
iROI = 3;
DATA.Name = char({'V Stim - A vs. V att';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,4); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(4,2,6)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,6)
iCond = 2;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0.5 0.5 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,8); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_V' ANTsSufix MedianSufix '.tif']), '-dtiff')
 
 
fclose (fid);

 
%% Plot Attention effect for V stim
DATA.Scatter = Scatter;
DATA.WithSubj = 1;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

close all

SavedTxt = fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_AV' ANTsSufix MedianSufix '.csv']);

fid = fopen (SavedTxt, 'w');
fprintf (fid, 'BOLD profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

figure('position', FigDim, 'name', 'MainEffect_Att_Stim_AV', 'Color', [1 1 1], 'visible', Visible)

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

DATA.MVPA = 0;
DATA.MIN = -1.2;
DATA.MAX = 2;
DATA.InsetLim = [1.2 0.4];
DATA.YLabelInset = 1;

DATA.OneSideTTest = {'right' 'left'};


% A1
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,1)
iCond = 3;
iROI = 1;
DATA.Name = char({'AV Stim - A vs. V att';'';'A1'});
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Color =  [1 0 0];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(BOLD_SubjectsBetas(:,1,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,2,:,2)),...
    squeeze(BOLD_SubjectsBetas(:,3,:,1)),...
    squeeze(BOLD_SubjectsBetas(:,4,:,1)));

PlotProfileAndBetas(DATA)

ax = subplot(4,2,3); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,5)
iCond = 3;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,2)
iCond = 3;
iROI = 3;
DATA.Name = char({'AV Stim - A vs. V att';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,4); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% V2-3
subplot(4,2,6)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,6)
iCond = 3;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0.5 0.5 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,8); 
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper',['MainEffect_Att_Stim_AV' ANTsSufix MedianSufix '.tif']), '-dtiff')
 
 
 
fclose (fid);
 