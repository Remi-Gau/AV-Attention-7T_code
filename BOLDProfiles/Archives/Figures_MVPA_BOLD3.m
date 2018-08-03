clc; clear; close all;

StartDirectory = fullfile(pwd, '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %     %'06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %     %'14';...
    '15';...
    '16'
    ];

% ROIs= {...
%     'A1_surf'
%     'PT_BT'
%     'V1_surf'
%     'V2-3_surf'
%     };

ROIs= {...
    'V1_A_Deact';...
    'V1_A_Act';...
    'V1_V_AV_Deact';...
    'V1_V_AV_Act';...
    'V2-3_A_Deact';...
    'V2-3_A_Act';...
    'V2-3_V_AV_Deact';...
    'V2-3_V_AV_Act';...
    
    };



% Holms Thresholds
HolmsThresholds = .05*ones(1,length(ROIs))./(length(ROIs)+1-(1:length(ROIs)));

FigDim = [100 100 1500 1000];
Visible = 'on';
Transparent = 1;
WithSubj = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);
NbLayers = 6;

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [DesMat' ones(NbLayers,1)];

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));


%% Get data for BOLD

% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 8 Interaction AV-A under A att vs. V under V att
% 9 Interaction AV-V under A att vs. V under V att
% 10 Interaction between A attention vs. V attention and A vs. V
% 11 A attention vs V attention - A Stim
% 12 A attention vs V attention - AV Stim
% 13 A attention vs V attention - V Stim
load(fullfile(FigureFolder, strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')
% load(fullfile(FigureFolder, strcat('Data2_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A against baseline & A against baseline
Target=1;
for iCond = 1:2
    for iROI=1:length(AllSubjects_Data)
        
        AllSubjects_Data(iROI).name
        AllSubjects_Data(iROI).VoxelCount
        AllSubjects_Data(iROI).ROI_Coverage(:,1)./AllSubjects_Data(iROI).ROI_Coverage(:,1)
        
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

% Interactions
Target=8;
for iCond = 1:3
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Interaction.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Interaction.DATA(:,iCond,iSubj);
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
% load(fullfile(FigureFolder, strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A_A vs A_V   &   V_A vs V_V
Target=3;
for iCond = 9:10
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Differential.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% AV-A and AV-V
Target=5;
for iCond = 7:8
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).BiVSUni.DATA(:,iCond,iSubj);
        end
    end
    Target = Target+1;
end

% Main effect of Att
Target=7;
for iCond = 2
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Differential.MainEffect.DATA(1:NbLayers,iCond,iSubj);
        end
    end
end


% Main effect of Att
Target=11;
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



%%
% Include =[];
% AllSubjects_Data(2).name='STGpost';
% for iROI=1:numel(ROIs)
%     if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
%         temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
%         Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
%     end
% end
% Include([4 11],:) = [];
% clear AllSubjects_Data temp iROI

Include = repmat(logical(ones(size(SubjectList,1),1)),[1,numel(ROIs)]);

% Holms Thresholds
for iCond = 1:size(All_P,3)
    for iP=1:size(All_P,2)
        tmp = All_P(:,iP,iCond);
        [~,I] = ismember(tmp,sort(tmp,'ascend'));
        BOLD_Thresholds(:,iP,iCond) = HolmsThresholds(I); %#ok<*SAGROW>
    end
end
clear All_P I iP iCond






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
% 11 A attention vs V attention - A Stim
% 12 A attention vs V attention - AV Stim
% 13 A attention vs V attention - V Stim

%% Plot
% Legends1 = {'', 'Constant', '', '', '', '', '', '', 'Linear', '', '', '', '', '', '', 'Quadratic',};
% Legends2 = {'mean', '(','SEM',')', 'p value', 'effect size', ''};
%
%
% SavedTxt = fullfile(StartDirectory, 'Figures', 'Deactivations_NoHolms_Quad.csv');
% fid1 = fopen (SavedTxt, 'w');
%
% fprintf (fid1, 'BOLD profile\n');
% fprintf (fid1, 'ROI');
% for i=1:length(Legends1)
%     fprintf (1, '%s,', Legends1{i});
% end
% fprintf (fid1, '\n,');
% for j=1:3
%     for i=1:length(Legends2)
%         fprintf (fid1, '%s,', Legends2{i});
%     end
% end
% fprintf (fid1, '\n');


% SavedTxt = fullfile(StartDirectory, 'Figures', 'Activations_NoHolms_Quad.csv');
% fid2 = fopen (SavedTxt, 'w');
%
% fprintf (fid2, 'BOLD profile\n');
% fprintf (fid2, 'ROI');
% for i=1:length(Legends1)
%     fprintf (fid2, '%s,', Legends1{i});
% end
% fprintf (fid2, '\n,');
% for j=1:3
%     for i=1:length(Legends2)
%         fprintf (fid2, '%s,', Legends2{i});
%     end
% end
% fprintf (fid2, '\n');



for iROI=1:numel(ROIs)
    
    % Plot BOLD
    for iCond = [1:2 5:7 11:13] %size(BOLD_SubjectsData,2)
        
        switch iCond
            case 1
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A against baseline irrespective of att';'(A_A + A_V)/2'});
            case 2
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'V against baseline irrespective of att';'(V_A + V_V)/2'});
            case 3
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A under A att vs A under V att';'A_A - A_V'});
            case 4
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'V under A att vs V under V att';'V_A - V_V'});
            case 5
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'AV-A irrespective of att';'( (AV-A)_A + (AV-A)_V )/2'});
            case 6
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'AV-V irrespective of att';'( (AV-V)_A + (AV-V)_V )/2'});
            case 7
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A att vs V att (irrespective of stim)';'( (A_A - A_V) + (V_A - V_V) + (AV_A - AV_V) )/3'});
            case 8
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'Interaction between AV vs A for A att vs V att';'(AV-A)_A - (AV-A)_V'});
            case 9
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'Interaction between AV vs V for A att vs V att';'(AV-V)_A - (AV-V)_V '});
            case 10
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'Interaction between A att vs V att and A vs V';'(A-V)_A - (A-V)_V'});
            case 11
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A att vs V att - A Stim';'(A_A - A_V)'});
            case 12
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A att vs V att - V Stim';'(V_A - V_V)'});
            case 13
                DATA.Name = char({strrep(ROIs{iROI},'_','-');'A att vs V att - AV Stim';'(AV_A - AV_V)'});
        end
        
        
        figure('position', FigDim, 'name', ROIs{iROI}, 'Color', [1 1 1], 'visible', Visible)
        
        set(gca,'units','centimeters')
        pos = get(gca,'Position');
        ti = get(gca,'TightInset');
        
        set(gcf, 'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        
        DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
        DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
        
        DATA.WithSubj = WithSubj;
        
        DATA.WithPerm = 0;
        DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
        DATA.Thresholds = ones(size(DATA.Thresholds))*0.05;
        
        DATA.XY= [0.65 0.65 .055 .15];
        DATA.Scatter = Scatter;
        
        DATA.MVPA = 0;
        
        DATA.YLabel = 'BOLD profile\n[a u]';
        
        PlotProfileAndBetasBasic(DATA)
        
        
        
        print(gcf, fullfile(StartDirectory, 'Figures', ...
            [ROIs{iROI} '_' num2str(iCond) '_' num2str(NbLayers) 'Layers_NoHolms_NoQuad.pdf']), '-dpdf')
        
        
        
        %          if (iCond==2 && iROI<3) || (iCond==1 && iROI>2)
        %              DATA.OneSideTTest = {'left' 'right' ''};
        %              fid = fid1;
        %          else
        %              DATA.OneSideTTest = {'right' 'left' ''};
        %              fid = fid2;
        %          end
        %
        %         fprintf (fid, '%s,', ROIs{iROI});
        %         for i=1:size(DATA.Betas,2)
        %             fprintf (fid, '%f,',nanmean(DATA.Betas(:,i)));
        %             fprintf (fid, '(,');
        %             fprintf (fid, '%f,',nansem(DATA.Betas(:,i)));
        %             fprintf (fid, '),');
        %             if ~isempty(DATA.OneSideTTest{i})
        %                 [~,P(i)] = ttest(DATA.Betas(:,i), 0, 'alpha', 0.05, 'tail', DATA.OneSideTTest{i});
        %             else
        %                 [~,P(i)] = ttest(DATA.Betas(:,i), 0, 'alpha', 0.05);
        %             end
        %             fprintf (fid, '%f,',P(i));
        %             fprintf (fid, '%f,',abs(nanmean(DATA.Betas(:,i))/nanstd(DATA.Betas(:,i))));
        %             fprintf (fid, ',');
        %         end
        %         fprintf (fid, '\n');
        
        clear DATA
        
    end
    
    close all
    
end

% fclose (fid1);
% fclose (fid2);





