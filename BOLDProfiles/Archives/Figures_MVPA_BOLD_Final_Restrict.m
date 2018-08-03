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
    %     'A1'
    %     'TE'
    %     'TE_surf'; ...
    'A1_surf'
    %     'STGpost'
    %     'PT_BT';...
    %     'pSTG_surf'; ...
    'PT_BT'
    'V1_surf'
    %     'V2-3'
    'V2-3_surf'
    };

Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');


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
load(fullfile(FigureFolder, strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A against baseline & A against baseline
Target=1;
for iCond = 1:2
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).MainEffectsRestrict.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% Interactions
Target=8;
for iCond = 1:3
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).InteractionRestrict.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end


% 3 A under A att vs. V under V att
% 4 V under A att vs. V under V att
% 5 AV-A irrespective of attention
% 6 AV-V irrespective of attention
% 7 A attention vs. V attention (pooled over all stimulation conditions)
load(fullfile(FigureFolder, strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A_A vs A_V   &   V_A vs V_V
Target=3;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).DifferentialRestrict.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% AV-A and AV-V
Target=5;
for iCond = 7:8
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).BiVSUniRestrict.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% Main effect of Att
Target=7;
for iCond = 2
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).DifferentialRestrict.MainEffect.DATA(1:NbLayers,iCond,iSubj);
        end
    end
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



%% Get Data for MVPA
cd(FigureFolder)

opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;
opt.scaling.idpdt = 1;


load(fullfile(FigureFolder, strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')


% Include =[];
% % AllSubjects_Data(2).name='STGpost';
% for iROI=1:numel(ROIs)
%     if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
%         temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
%         Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
%     end
% end
% Include([4 11],:) = [];
% clear AllSubjects_Data temp iROI

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
            %             SaveSufix = [SaveSufix '_FWHM_6Slab_Layers_' num2str(NbLayers) '.mat'];
            
            Save_vol = [Save_vol SaveSufix];
            
            try
                load(fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID],  'Transfer', 'SVM', ...
                    Save_vol), 'Class_Acc', 'Results');
                MVPA_SubjectsData(SubjInd,iSVM,iROI,1:(NbLayers+1)) = squeeze(Class_Acc.TotAcc(end,:,:,:));
                
                for iLayer = 2:NbLayers+1
                    Label = cat(2,Results(end).session(end).rand.perm.SubSamp{1}.CV(:,iLayer).label);
                    Pred = cat(2,Results(end).session(end).rand.perm.SubSamp{1}.CV(:,iLayer).pred);
                    Acc(iLayer-1,:) = mean(Pred==Label,2)';
                end
                
                X=repmat(DesMat,size(Acc,2),1);
                Acc=flipud(Acc(:)-.5);
                [B,~,~] = glmfit(X, Acc, 'normal', 'constant', 'off');
                SubjectsBetas(SubjInd,iSVM,iROI,1:size(X,2)) = B;
                
            catch
                warning('file %s missing', fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID],  'Transfer', 'SVM', ...
                    Save_vol))
                MVPA_SubjectsData(SubjInd,iSVM,iROI,1:(NbLayers+1)) = NaN;
                SubjectsBetas(SubjInd,iSVM,iROI,:) = NaN;
            end
            
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

% Holms Thresholds
for iSVM = 1:numel(Analysis)
    for iP=1:size(All_P,2)
        tmp = All_P(:,iP,iSVM);
        [~,I] = ismember(tmp,sort(tmp,'ascend'));
        MVPA_Thresholds(:,iP,iSVM) = HolmsThresholds(I);
    end
end
clear All_P






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


Legends1 = {'', 'Constant', '','','', '', '', '', 'Linear'};
Legends2 = {'ROI', 'mean', '(','SEM',')', 'p value', 'effect size', '', 'mean', '(','SEM',')', 'p value', 'effect size'};


%% Plot deactivations

% SavedTxt = fullfile(StartDirectory, 'Paper','Deactivations.csv');
% fid = fopen (SavedTxt, 'w');

DATA.WithSubj = 1;

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

% fprintf (fid, 'BOLD profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');

DATA.MVPA = 0;
DATA.YLabelInset = 1;
DATA.InsetLim = [0.6 0.18];
if DATA.WithSubj
    DATA.MIN = -0.6;
    DATA.MAX = 0.45;
else
    DATA.MIN = -0.3;
    DATA.MAX = 0.1;
end

DATA.OneSideTTest = {'left' 'right'};


% A1
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,1)
iCond = 2;
iROI = 1;
DATA.Name = char({'V vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,5)
iCond = 2;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



DATA.InsetLim = [1.6 0.35];
if DATA.WithSubj
    DATA.MIN = -2.2;
    DATA.MAX = 0.1;
else
    DATA.MIN = -1.5;
    DATA.MAX = 0.1;
end

% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,2)
iCond = 1;
iROI = 3;
DATA.Name = char({'A vs. Baseline';'';'V1'});
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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
    'Fig2_Deactivations_Restrict.tif'), '-dtiff')

% plot2svg(fullfile(StartDirectory, 'Paper', ...
%      'Fig2_Deactivations.svg'),gcf, 'svg')

% fclose (fid);

% saveas(gcf, fullfile(StartDirectory, 'Paper', ...
%       'Fig2_Deactivations.eps'),'epsc2')
%
% % 'svg'



%% Plot Activations
clear DATA

% SavedTxt = fullfile(StartDirectory, 'Paper','Activations.csv');
% fid = fopen (SavedTxt, 'w');

DATA.WithSubj = 1;

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

% fprintf (fid, 'BOLD profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');

DATA.MVPA = 0;
DATA.YLabelInset = 1;

DATA.InsetLim = [2.8 1.1];
if DATA.WithSubj
    DATA.MIN = -0.5;
    DATA.MAX = 5.5;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end


DATA.OneSideTTest = {'right' 'left'};


% TE
subplot(4,2,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,1)
iCond = 1;
iROI = 1;
DATA.Name = char({'A vs. Baseline';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)



% PT
subplot(4,2,5)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,5)
iCond = 1;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,2,7);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)

DATA.InsetLim = [2.8 1.1];
if DATA.WithSubj
    DATA.MIN = -0.1;
    DATA.MAX = 4.4;
else
    DATA.MIN = -0.1;
    DATA.MAX = 3.1;
end

% V1
subplot(4,2,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,2,2)
iCond = 2;
iROI = 3;
DATA.Name = char({'V vs. Baseline';'';'V1'});
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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


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

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
    'Fig1_Activations_Restrict.tif'), '-dtiff')

% fclose (fid);

%% Plot Cross Modal Influence
clear DATA

% SavedTxt = fullfile(StartDirectory, 'Paper','CrossModal.csv');
% fid = fopen (SavedTxt, 'w');

DATA.WithSubj = 1;

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
DATA.InsetLim = [1 0.3];
if DATA.WithSubj
    DATA.MIN = -1.3;
    DATA.MAX = 2;
else
    DATA.MIN = -0.25;
    DATA.MAX = 0.7;
end

% fprintf (fid, 'BOLD profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');


% TE
subplot(4,4,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,1)
iCond = 5;
iROI = 1;
DATA.Name = char({'AV vs. A';'';'A1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0 0];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(BOLD_SubjectsBetas(:,1,:,5)),...
    squeeze(BOLD_SubjectsBetas(:,2,:,5)),...
    squeeze(BOLD_SubjectsBetas(:,3,:,6)),...
    squeeze(BOLD_SubjectsBetas(:,4,:,6)));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,5);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 1;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


% PT
subplot(4,4,9)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,9)
iCond = 5;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,13);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)

DATA.InsetLim = [1 0.3];
if DATA.WithSubj
    DATA.MIN = -2;
    DATA.MAX = 1;
else
        DATA.MIN = -0.51;
    DATA.MAX = 0.11;
end

% V1
subplot(4,4,3)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,3)
iCond = 6;
iROI = 3;
DATA.Name = char({'AV vs. V';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,7);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


% V2-3
subplot(4,4,11)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,11)
iCond = 6;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0.5 0.5 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,15);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)



% Plot MVPA
DATA.MVPA = 1;
DATA.InsetLim = [0.3 0.075];
if DATA.WithSubj
    DATA.MIN = 0.3;
    DATA.MAX = 0.95;
else
    DATA.MIN = 0.49;
    DATA.MAX = 0.72;
end

DATA.YLabel = 'Accuracy';

% fprintf (fid, '\n\n');
% fprintf (fid, 'MVPA profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');

% TE
subplot(4,4,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,2)
iSVM = 1;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [1 0 0];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(MVPA_SubjectsBetas(:,1,:,1)),...
    squeeze(MVPA_SubjectsBetas(:,2,:,1)),...
    squeeze(MVPA_SubjectsBetas(:,3,:,2)),...
    squeeze(MVPA_SubjectsBetas(:,4,:,2)));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,6);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)

% PT
subplot(4,4,10)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,10)
iSVM = 1;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [1 0.5 0.5];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,14);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)

DATA.MVPA = 1;
DATA.InsetLim = [0.3 0.05];
    
if DATA.WithSubj
    DATA.MIN = 0.3;
    DATA.MAX = 0.83;
else
    DATA.MIN = 0.44;
    DATA.MAX = 0.64;
end


% V1
subplot(4,4,4)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,4)
iSVM = 2;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [0 0 1];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,8);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


% V2-3
subplot(4,4,12)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,12)
iSVM = 2;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [0.5 0.5 1];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,16);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
    'Fig3_CrossModalInfluence_Restrict.tif'),  '-dtiff')

% fclose (fid);


%% Plot attention effects
clear DATA

% SavedTxt = fullfile(StartDirectory, 'Paper','Attention.csv');
% fid = fopen (SavedTxt, 'w');

DATA.WithSubj = 1;

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
DATA.InsetLim = [0.8 0.2];
if DATA.WithSubj
    DATA.MIN = -0.6;
    DATA.MAX = 1.1;
else
    DATA.MIN = -0.1;
    DATA.MAX = 0.6;
end

DATA.OneSideTTest = {'right' 'left'};

% fprintf (fid, 'BOLD profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');

% TE
subplot(4,4,1)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,1)
iCond = 7;
iROI = 1;
DATA.Name = char({'A vs. V att';'';'A1'});
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0 0];
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(BOLD_SubjectsBetas(:,1,:,7)),...
    squeeze(BOLD_SubjectsBetas(:,2,:,7)),...
    squeeze(BOLD_SubjectsBetas(:,3,:,7)),...
    squeeze(BOLD_SubjectsBetas(:,4,:,7)));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,5);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 1;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% PT
subplot(4,4,9)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,9)
iCond = 7;
iROI = 2;
DATA.Name = 'PT';
DATA.Data = -1*fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [1 0.5 0.5];
DATA.Betas = -1*squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,13);
axis('off')
DATA.ax = ax.Position;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


if DATA.WithSubj
    DATA.MIN = -0.6;
    DATA.MAX = 1.1;
else
    DATA.MIN = -0.1;
    DATA.MAX = 0.4;
end

% V1
subplot(4,4,3)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,3)
iCond = 7;
iROI = 3;
DATA.Name = char({'V vs. A att';'';'V1'});
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0 0 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,7);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% V2-3
subplot(4,4,11)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,11)
iCond = 7;
iROI = 4;
DATA.Name = 'V2-3';
DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
DATA.Color =  [0.5 0.5 1];
DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,15);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA, DATA.OneSideTTest)


% Plot MVPA
DATA = rmfield(DATA, 'OneSideTTest');

DATA.MVPA = 1;
DATA.YLabel = 'Accuracy';

DATA.InsetLim = [0.14 0.05];
if DATA.WithSubj
    DATA.MIN = 0.35;
    DATA.MAX = 0.75;
else
    DATA.MIN = 0.45;
    DATA.MAX = 0.62;
end

% fprintf (fid, '\n\n');
% fprintf (fid, 'MVPA profile\n');
% for i=1:length(Legends1)
%     fprintf (fid, '%s,', Legends1{i});
% end
% fprintf (fid, '\n');
% for i=1:length(Legends2)
%     fprintf (fid, '%s,', Legends2{i});
% end
% fprintf (fid, '\n');

% TE
subplot(4,4,2)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,2)
iSVM = 3;
iROI = 1;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [1 0 0];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));

DATA.OtherBetas = cat(3,...
    squeeze(MVPA_SubjectsBetas(:,1,:,3)),...
    squeeze(MVPA_SubjectsBetas(:,2,:,3)),...
    squeeze(MVPA_SubjectsBetas(:,3,:,3)),...
    squeeze(MVPA_SubjectsBetas(:,4,:,3)));

PlotProfileAndBetas(DATA)

ax = subplot(4,4,6);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


% PT
subplot(4,4,10)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,10)
iSVM = 3;
iROI = 2;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [1 0.5 0.5];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,14);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)

DATA.InsetLim = [0.28 0.06];
if DATA.WithSubj
    DATA.MIN = 0.3;
    DATA.MAX = 0.81;
else
    DATA.MIN = 0.49;
    DATA.MAX = 0.67;
end

% V1
subplot(4,4,4)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,4)
iSVM = 3;
iROI = 3;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [0 0 1];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,8);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


% V2-3
subplot(4,4,12)
PlotRectangle(NbLayers,Fontsize)
subplot(4,4,12)
iSVM = 3;
iROI = 4;
DATA.Name = '';
DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
DATA.Color =  [0.5 0.5 1];
DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
DATA.Thresholds = 0.05*ones(size(DATA.Thresholds));
PlotProfileAndBetas(DATA)

ax = subplot(4,4,16);
axis('off')
DATA.ax = ax.Position;
DATA.YLabelInset = 0;
PlotInset(DATA)

% Print2Table(fid, ROIs, iROI, DATA)


cd(fullfile(StartDirectory, 'Paper'))
print(gcf, fullfile(StartDirectory, 'Paper', ...
    'Fig4_Attention_Restrict.tif'), '-dtiff')

% fclose (fid);