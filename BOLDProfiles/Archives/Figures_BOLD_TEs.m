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
    'TE1.0_surf';...
    'TE1.1_surf';...
    'TE1.2_surf';...
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
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = [DesMat' ones(NbLayers,1)];
DesMat = spm_orth(DesMat);

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));


%% Get data for BOLD

% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
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
        end
    end
    Target = Target+1;
end



%% Get data for BOLD 2

% 3 AV-A irrespective of attention
% 4 AV-V irrespective of attention
% 5 A attention vs. V attention (pooled over all stimulation conditions)
load(fullfile(FigureFolder, strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')
% load(fullfile(FigureFolder, strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')


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
            BOLD_SubjectsData(iSubj,Target,iROI,1:NbLayers) = AllSubjects_Data(iROI).Differential.MainEffect.DATA(1:NbLayers,iCond,iSubj);
        end
    end
end



%%
Include = repmat(logical(ones(size(SubjectList,1),1)),[1,numel(ROIs)]);


% 1 A against baseline irrespective of attention
% 2 V against baseline irrespective of attention
% 3 AV-A irrespective of attention
% 4 AV-V irrespective of attention
% 5 A attention vs V attention (pooled over all stimulation conditions)

%% Plot

for iROI=1:numel(ROIs)
    
    % Plot BOLD
    for iCond = 1:size(BOLD_SubjectsData,2)
        
        switch iCond
            case 1
                DATA.Name = char({['PSC - ' strrep(ROIs{iROI},'_', ' ')];'A against baseline irrespective of att';'(A_A + A_V)/2'});
            case 2
                DATA.Name = char({['PSC - ' strrep(ROIs{iROI},'_', ' ')];'V against baseline irrespective of att';'(V_A + V_V)/2'});
            case 3
                DATA.Name = char({['PSC - ' strrep(ROIs{iROI},'_', ' ')];'AV-A irrespective of att';'( (AV-A)_A + (AV-A)_V )/2'});
            case 4
                DATA.Name = char({['PSC - ' strrep(ROIs{iROI},'_', ' ')];'AV-V irrespective of att';'( (AV-V)_A + (AV-V)_V )/2'});
            case 5
                DATA.Name = char({['PSC - ' strrep(ROIs{iROI},'_', ' ')];'A att vs V att (irrespective of stim)';'( (A_A - A_V) + (V_A - V_V) + (AV_A - AV_V) )/3'});
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
        
        DATA.XY= [0.65 0.65 .055 .15];
        DATA.Scatter = Scatter;
        
        DATA.MVPA = 0;

        DATA.YLabel = 'Percent signal change';

        PlotProfileAndBetasBasic(DATA)
        clear DATA

        print(gcf, fullfile(StartDirectory, 'Figures', ...
            [ROIs{iROI} '_' num2str(iCond) '_' num2str(NbLayers) 'Layers_NoQuad.pdf']), '-dpdf')

        
    end
    
    close all
    
end





