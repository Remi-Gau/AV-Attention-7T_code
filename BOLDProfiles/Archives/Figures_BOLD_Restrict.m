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

ROIs= {...
    'A1_surf'
    'PT_BT'
    'V1_surf'
    'V2-3_surf'
    };

FigDim = [100 100 1500 1000];
Visible = 'on';
Transparent = 1;
WithSubj = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);
NbLayers = 6;

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));


%% Get data for BOLD

% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 8 Interaction AV-A under A att vs. V under V att
% 9 Interaction AV-V under A att vs. V under V att
% 10 Interaction between A attention vs. V attention and A vs. V

load(fullfile(FigureFolder, strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')
% load(fullfile(FigureFolder, strcat('Data2_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

% A against baseline & A against baseline
Target=1;
for iCond = 1:2
    for iROI=1:length(AllSubjects_Data)
        Include = find(AllSubjects_Data(iROI).Include);
        tmp = squeeze(AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA(:,iCond,Include))';
        [~,P] = ttest(tmp);
        All_P(iROI,:,Target) = P;
        BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),Target) = tmp;
        
        for iSubj = 1:length(Include)
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).MainEffects.DATA(:,iCond,iSubj);
            BOLD_Restrict_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).MainEffectsRestrict.DATA(:,iCond,iSubj);
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
            BOLD_Restrict_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).InteractionRestrict.DATA(:,iCond,iSubj);
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
            BOLD_Restrict_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).DifferentialRestrict.DATA(:,iCond,iSubj);
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
            BOLD_Restrict_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).BiVSUniRestrict.DATA(:,iCond,iSubj);
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
            BOLD_Restrict_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).DifferentialRestrict.DATA(:,iCond,iSubj);
        end
    end
end



%%
% Include = true(size(SubjectList,1));

Include = false(size(SubjectList,1));
Include(7) = true;

Include = repmat(Include,[1,numel(ROIs)]);

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

%% Plot

for iROI=1:numel(ROIs)
    
    % Plot BOLD
    for iCond = [1 2 5 6 7] %1:size(BOLD_SubjectsData,2)
        
        switch iCond
            case 1
                DATA.Name = char({ROIs{iROI};'A against baseline irrespective of att';'(A_A + A_V)/2'});
            case 2
                DATA.Name = char({ROIs{iROI};'V against baseline irrespective of att';'(V_A + V_V)/2'});
            case 3
                DATA.Name = char({ROIs{iROI};'A under A att vs A under V att';'A_A - A_V'});
            case 4
                DATA.Name = char({ROIs{iROI};'V under A att vs V under V att';'V_A - V_V'});
            case 5
                DATA.Name = char({ROIs{iROI};'AV-A irrespective of att';'( (AV-A)_A + (AV-A)_V )/2'});
            case 6
                DATA.Name = char({ROIs{iROI};'AV-V irrespective of att';'( (AV-V)_A + (AV-V)_V )/2'});
            case 7
                DATA.Name = char({ROIs{iROI};'A att vs V att (irrespective of stim)';'( (A_A - A_V) + (V_A - V_V) + (AV_A - AV_V) )/3'});
            case 8
                DATA.Name = char({ROIs{iROI};'Interaction between AV vs A for A att vs V att';'(AV-A)_A - (AV-A)_V'});
            case 9
                DATA.Name = char({ROIs{iROI};'Interaction between AV vs V for A att vs V att';'(AV-V)_A - (AV-V)_V '});
            case 10
                DATA.Name = char({ROIs{iROI};'Interaction between A att vs V att and A vs V';'(A-V)_A - (A-V)_V'});
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
        DATA.DataRestrict = fliplr(squeeze(BOLD_Restrict_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));

        DATA.WithSubj = WithSubj;

        DATA.MVPA = 0;

        DATA.YLabel = 'BOLD profile\n[a u]';

        PlotProfileAndBetasBasicRestrict(DATA)
        clear DATA
        
%         print(gcf, fullfile(StartDirectory, 'Figures', ...
%             [ROIs{iROI} '_' num2str(iCond) '_' num2str(NbLayers) 'Layers_NoHolms_Quad.pdf']), '-dpdf')

        
    end
    
    close all
    
end





