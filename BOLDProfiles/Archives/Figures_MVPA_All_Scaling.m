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


% Analysis(1) = struct('name', 'A Stim VS V Stim');
% Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'AV VS A+V');
% Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');

ROIs= {...
%     'TE'
%     'STGpost'
    'V1'
%     'V2-3'
    %'V3'
    %     'V4'
    %     'V5'
    %     'V1-2-3'
    };

% Holms Thresholds
HolmsThresholds = .05*ones(1,length(ROIs))./(length(ROIs)+1-(1:length(ROIs)));

FigDim = [100 100 1500 1000];

Visible = 'on';

Transparent = 1;

Fontsize = 12;

Scatter = linspace(0,.4,11);

NbLayers = 6;

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));


opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;
opt.scaling.idpdt = 1;
WithSubj = 1;


DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [DesMat' (DesMat.^2)' ones(NbLayers,1)];
% DesMat = [DesMat' ones(NbLayers,1)];



%% Get Data for MVPA
cd(FigureFolder)



load(fullfile(FigureFolder, strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')


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
            
            SaveSufix = [SaveSufix '_FWHM_6Slab_Layers_' num2str(NbLayers) '.mat'];
            
            Save_vol = [Save_vol SaveSufix];
            
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


%% Plot
for iROI=1:numel(ROIs)
    
    for iSVM = 1:numel(Analysis)
        
        figure('position', FigDim, 'name', ROIs{iROI}, 'Color', [1 1 1], 'visible', Visible)
        
        set(gca,'units','centimeters')
        pos = get(gca,'Position');
        ti = get(gca,'TightInset');
        
        set(gcf, 'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        
        DATA.Data = fliplr(squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
        DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
        
        DATA.WithSubj = WithSubj;
        
        DATA.WithPerm = 0;
        DATA.Thresholds = squeeze(MVPA_Thresholds(iROI,:,iSVM));
        DATA.Thresholds = ones(size(DATA.Thresholds))*0.05;
        
        DATA.XY= [0.65 0.65 .055 .15];
        DATA.Scatter = Scatter;
        
        DATA.MVPA = 1;
        
        DATA.YLabel = 'Decoding profile\n[Percent correct]';
        
        DATA.Name = char({ROIs{iROI};['MVPA: ' Analysis(iSVM).name]});
        
        PlotProfileAndBetasBasic(DATA)
        clear DATA

%         print(gcf, fullfile(StartDirectory, 'Figures', ...
%             [ROIs{iROI} '_MVPA_' num2str(iSVM) '_' num2str(NbLayers) 'Layers_NoHolms.pdf']), '-dpdf')
%         
                print(gcf, fullfile(StartDirectory, 'Figures', ...
            [ROIs{iROI} '_MVPA_' num2str(iSVM) '_' num2str(NbLayers) 'Layers_NoHolms_Quad.pdf']), '-dpdf')
        
%                 print(gcf, fullfile(StartDirectory, 'Figures', ...
%             [ROIs{iROI} '_MVPA_' num2str(iSVM) '_' num2str(NbLayers) 'Layers_NoHolms.tif']), '-dtiff')

    end
    
    close all
    
end
