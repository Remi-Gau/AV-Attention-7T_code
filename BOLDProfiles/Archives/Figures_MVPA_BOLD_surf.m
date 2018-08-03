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

% Get all possible permutations
for iSubj=1:size(SubjectList,1)
    sets{iSubj} = [0 1];
end
[a,b,c,d,e,f,g,h,i,j,k] = ndgrid(sets{:});
PermList = logical([a(:),b(:),c(:),d(:),e(:),f(:),g(:),h(:),i(:),j(:),k(:)]);

% [a,b,c,d,e,f,g,h,i,j] = ndgrid(sets{:});
% PermList = logical([a(:),b(:),c(:),d(:),e(:),f(:),g(:),h(:),i(:),j(:)]);

NbPerm=size(PermList,1);
clear sets a b c d e f g h i j k iSubj

% Analysis(1) = struct('name', 'A Stim VS V Stim');
Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'AV VS A+V');
Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');

ROIs= {...
    'TE'
    'STGpost'
    'V1'
    'V2-3'
    %'V3'
    %     'V4'
    %     'V5'
    %     'V1-2-3'
    };

% Holms Thresholds
HolmsThresholds = .05*ones(1,length(ROIs))./(length(ROIs)+1-(1:length(ROIs)));

FigDim = [100 100 1500 1000];

Visible = 'on';

n=2;
m=4;

Transparent = 1;

Fontsize = 12;

Scatter = linspace(0,.4,11);

for NbLayers = [8]
    
    FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers-2), '_layers'));
    
    for ImgNorm= 0
        
        opt.scaling.img.eucledian = 0;
        opt.scaling.img.zscore = ImgNorm;
        opt.scaling.feat.mean = 1;
        opt.scaling.feat.range = 0;
        opt.scaling.feat.sessmean = 0;
        opt.scaling.idpdt = 1;
        
        for WithQuad= 0
            
            DesMat = (1:(NbLayers-2))-mean(1:(NbLayers-2));
            if WithQuad
                DesMat = [DesMat' (DesMat.^2)' ones((NbLayers-2),1)];
            else
                DesMat = [DesMat' ones((NbLayers-2),1)];
            end
            
            
            
            %% Get data for BOLD
            if WithQuad
                load(fullfile(FigureFolder, strcat('Data_Block_Surf_', num2str(NbLayers-2), '_Layers', '.mat')), 'AllSubjects_Data') %#ok<*UNRCH>
            else
                load(fullfile(FigureFolder, strcat('Data_Block_Surf_NoQuadGLM_', num2str(NbLayers-2), '_Layers', '.mat')), 'AllSubjects_Data')
            end
            
            % Compile Betas and p values for AV-A and AV-V
            for iCond = 1:2
                for iROI=1:length(AllSubjects_Data)
                    Include = find(AllSubjects_Data(iROI).Include);
                    tmp = squeeze(AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.DATA(:,iCond+6,Include))';
                    [~,P] = ttest(tmp);
                    All_P(iROI,:,iCond) = P;
                    BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),iCond) = tmp;
                    
                    for iSubj = 1:length(Include)
                        BOLD_SubjectsData(iSubj,iCond,iROI,1:(NbLayers-2)) = ...
                            AllSubjects_Data(iROI).BiVSUni.DATA(:,iCond+6,iSubj);
                    end
                end
            end
            clear P
            
            
            % Compile Betas and p values for main effect of Att and MSI
            for iCond = 1:2
                for iROI=1:length(AllSubjects_Data)
                    Include = find(AllSubjects_Data(iROI).Include);
                    tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA(:,iCond,Include))';
                    [~,P] = ttest(tmp);
                    All_P(iROI,:,iCond+2) = P;
                    BOLD_SubjectsBetas(:,iROI,1:size(tmp,2),iCond+2) = tmp;
                    
                    for iSubj = 1:length(Include)
                        BOLD_SubjectsData(iSubj,iCond+2,iROI,1:(NbLayers-2)) = ...
                            AllSubjects_Data(iROI).Differential.MainEffect.DATA(1:(NbLayers-2),iCond,iSubj);
                    end
                end
            end
            clear P
            
            
            
            % Holms Thresholds
            for iCond = 1:size(All_P,3)
                for iP=1:size(All_P,2)
                    tmp = All_P(:,iP,iCond);
                    [~,I] = ismember(tmp,sort(tmp,'ascend'));
                    BOLD_Thresholds(:,iP,iCond) = HolmsThresholds(I); %#ok<*SAGROW>
                end
            end
            clear All_P I iP iCond
            
            % Do the permutations
            NullDist = nan(NbPerm,size(BOLD_SubjectsBetas,3));
            for iCond = 1:size(BOLD_SubjectsBetas,4)
                for iPerm=1:NbPerm
                    tmp = BOLD_SubjectsBetas(:,:,:,iCond);
                    tmp(PermList(iPerm,:),:,:) = tmp(PermList(iPerm,:),:,:)*-1;
                    NullDist(iPerm,:) = squeeze(max(abs(mean(tmp)),[],2));
                end
                BOLD_AllNullDist(:,:,iCond) = sort(NullDist);
                clear NullDist
            end
            clear All_Betas iPerm tmp AllSubjects_Data
            
            
            %% Get Data for MVPA
            cd(FigureFolder)
            

            if WithQuad
                load(fullfile(FigureFolder, strcat('Data_Block_Surf_', num2str(NbLayers-2), '_Layers', '.mat')), 'AllSubjects_Data') %#ok<*UNRCH>
            else
                load(fullfile(FigureFolder, strcat('Data_Block_Surf_NoQuadGLM_', num2str(NbLayers-2), '_Layers', '.mat')), 'AllSubjects_Data')
            end
            
            
            Include =[];
            AllSubjects_Data(2).name='STGpost';
            for iROI=1:numel(ROIs)
                if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
                    temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
                    Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
                end
            end
            clear AllSubjects_Data temp iROI
            
            
            for SubjInd = 1:size(SubjectList,1)
                
                SubjID = SubjectList(SubjInd,:);
                
                %%
                for iSVM=1:numel(Analysis)
                    
                    for iROI=1:numel(ROIs)
                        
                        if strcmp(ROIs{iROI}, 'STGpost')
                            Save_vol = [...
                                'SVM_' Analysis(iSVM).name...
                                '_ROI_pSTG'];
                        else
                                                    Save_vol = [...
                            'SVM_' Analysis(iSVM).name...
                            '_ROI_' ROIs{iROI}];
                        end
                        
                        SaveSufix = '_results_surf_FixedC';
                        
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
                        
                        SaveSufix = [SaveSufix '_FWHM_6_Layers_' num2str(NbLayers) '.mat'];
                        
                        Save_vol = [Save_vol SaveSufix];
                        
                        load(fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID],  'Transfer', 'SVM', ...
                            Save_vol), 'Results');
                        
                        % get accuracy data
                        for iSes=1:numel(Results.session(end).rand.perm.CV)
                            tmp(iSes,:)=diag(Results.session(end).rand.perm.CV(iSes).layers.acc);
                        end
                        MVPA_SubjectsData(SubjInd,iSVM,iROI,1:(NbLayers-2)) = mean(tmp(:,2:end-1));
                        clear tmp
                        
                        % Get prediction for each block
                        for iLayer = 2:NbLayers-1
                            Label = Results.session(end).rand.perm.CV(iSes).layers.results{1}{iLayer}.label;
                            Pred = Results.session(end).rand.perm.CV(iSes).layers.results{1}{iLayer}.pred(:,iLayer);
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
            
            % Do the permutations
            NullDist = nan(NbPerm,size(MVPA_SubjectsBetas,3));
            for iSVM = 1:numel(Analysis)
                for iPerm=1:NbPerm
                    tmp = MVPA_SubjectsBetas(:,:,:,iSVM);
                    tmp(PermList(iPerm,:),:,:) = tmp(PermList(iPerm,:),:,:)*-1;
                    NullDist(iPerm,:) = squeeze(max(abs(mean(tmp)),[],2));
                end
                MVPA_AllNullDist(:,:,iSVM) = sort(NullDist);
            end
            
            
            %% Plot
            
            for WithPerm = 0
                
                for WithSubj=[1]
                    
                    for iROI=1:numel(ROIs)
                        
                        figure('position', FigDim, 'name', ROIs{iROI}, 'Color', [1 1 1], 'visible', Visible)
                        
                        set(gca,'units','centimeters')
                        pos = get(gca,'Position');
                        ti = get(gca,'TightInset');
                        
                        set(gcf, 'PaperUnits','centimeters');
                        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                        set(gcf, 'PaperPositionMode', 'manual');
                        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                        
                        title(ROIs{iROI})
                        
                        % Plot BOLD
                        iSubPlot=1;
                        
                        for iCond = 1:size(BOLD_SubjectsData,2)
                            
                            subplot(n,m,iSubPlot)
                            
                            switch iSubPlot
                                case 1
                                    XY = [0.17 0.625 .04 .07];
                                    DATA.Name = '([AV-A)]_A + [AV-A)]_V)/2';
                                case 2
                                    XY = [0.375 0.625 .04 .07];
                                    DATA.Name = '([AV-V)]_A + [AV-V)]_V)/2';
                                case 3
                                    XY = [0.58 0.625 .04 .07];
                                    DATA.Name = 'Main effect of MSI';
                                case 4
                                    XY = [0.778 0.625 .04 .07];
                                    DATA.Name = 'Main effect of attention';
                            end
                            % 'Main effect of MSI: $\frac{[AV-(A+V)]_{A}+[AV-(A+V)]_{V}}{2}$';...$\frac{[AV-(A+V)]_{A}+[AV-(A+V)]_{V}}{2}$ ([AV-(A+V)]_A + [AV-(A+V)]_V)/2
                            % 'Main effect of attention: $\frac{[Att_V-Att_A]_{V}+[Att_V-Att_A]_{A}+[Att_V-Att_A]_{AV}}{3}$'}; %$\frac{[Att_V-Att_A]_{V}+[Att_V-Att_A]_{A}+[Att_V-Att_A]_{AV}}{3}$' ([Att_V-Att_A]_V + [Att_V-Att_A]_A + [Att_V-Att_A]_{AV})/3'
                            
                            
                            DATA.Data = fliplr(squeeze(BOLD_SubjectsData(logical(Include(:,iROI)),iCond,iROI,:)));
                            DATA.Betas = squeeze(BOLD_SubjectsBetas(:,iROI,:,iCond));
                            
                            DATA.WithSubj = WithSubj;
                            
                            DATA.WithPerm = WithPerm;
                            if WithPerm
                                DATA.AllNullDist = BOLD_AllNullDist(:,:,iCond);
                            else
                                DATA.Thresholds = squeeze(BOLD_Thresholds(iROI,:,iCond));
                            end
                            
                            DATA.XY= XY;
                            DATA.Scatter = Scatter;
                            
                            DATA.MVPA = 0;
                            
                            if iSubPlot==1
                                DATA.YLabel = 'BOLD profile\n[a u]';
                            end
                            
                            PlotProfileAndBetas(DATA)
                            clear DATA
                            
                            iSubPlot = iSubPlot + 1;
                            
                        end
                        
                        
                        % Plot MVPA
                        iSubPlot=5;
                        
                        for iSVM = 1:numel(Analysis)
                            
                            subplot(n,m,iSubPlot)
                            
                            DATA.Name = Analysis(iSVM).name;
                            
                            DATA.Data = squeeze(MVPA_SubjectsData(logical(Include(:,iROI)),iSVM,iROI,1:(NbLayers-2)));
                            DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,:,iSVM));
                            
                            DATA.WithSubj = WithSubj;
                            
                            DATA.WithPerm = WithPerm;
                            if WithPerm
                                DATA.AllNullDist = MVPA_AllNullDist(:,:,iSVM);
                            else
                                DATA.Thresholds = squeeze(MVPA_Thresholds(iROI,:,iSVM));
                            end
                            
                            
                            switch iSubPlot
                                case 5
                                    XY = [0.17 0.15 .04 .07];
                                case 6
                                    XY = [0.375 0.15 .04 .07];
                                case 7
                                    XY = [0.58 0.15 .04 .07];
                            end
                            
                            DATA.XY= XY;
                            DATA.Scatter = Scatter;
                            
                            DATA.MVPA = 1;
                            
                            if iSubPlot==5
                                DATA.YLabel = 'Decoding profile\n[Percent correct]';
                            end
                            
                            PlotProfileAndBetas(DATA)
                            clear DATA
                            
                            iSubPlot = iSubPlot + 1;
                            
                        end
                        
                        
                        %%
                        SaveSufix = 'SURFACE\n';
                        
                        if opt.scaling.idpdt
                            SaveSufix = [SaveSufix 'Idpdt Norm\n']; %#ok<*AGROW>
                        end
                        if opt.scaling.img.zscore
                            SaveSufix = [SaveSufix 'IMG norm: ZScore\n'];
                        end
                        if opt.scaling.img.eucledian
                            SaveSufix = [SaveSufix 'IMG norm: Eucl\n'];
                        end
                        if opt.scaling.img.zscore==0 && opt.scaling.img.eucledian==0
                            SaveSufix = [SaveSufix 'IMG norm: None\n'];
                        end
                        if opt.scaling.feat.mean
                            SaveSufix = [SaveSufix 'FEAT norm: MeanCent\n'];
                        end
                        if opt.scaling.feat.range
                            SaveSufix = [SaveSufix 'FEAT norm: Range\n'];
                        end
                        if opt.scaling.feat.sessmean
                            SaveSufix = [SaveSufix 'FEAT norm: SessMeanCent\n'];
                        end
                        SaveSufix = [SaveSufix 'Smoothing: FWHM 6 mm\n'];
                        
                        subplot(n,m,8)
                        t=title(ROIs{iROI});
                        set(t,'fontsize', Fontsize+4)
                        axis off
                        t=text(0.1, 0.8, sprintf(SaveSufix));
                        set(t,'fontsize', Fontsize)
                        
                        Name = strrep([ROIs{iROI} '_' SaveSufix], '\n', '_');
                        Name = strrep(Name, ': ', '_');
                        Name = strrep(Name, ' ', '_');
                        
                        if WithQuad
                            Name = [Name 'QUAD_'];
                        end
                        if WithPerm
                            Name = [Name 'PERM_'];
                        end
                        if WithSubj
                            Name = [Name 'SUBJ_'];
                        end
                        
                        
                        print(gcf, fullfile(StartDirectory, 'Figures', ...
                            [Name num2str(NbLayers) 'Layers_surf.pdf']), '-dpdf')
                        
                        %                         print(gcf, fullfile(StartDirectory, 'Figures', 'MVPA', 'vol', ...
                        %                             [Name num2str(NbLayers) 'Layers.tif']), '-dtiff')
                        
                        close all
                        
                    end
                    
                    % clear AllNullDist AllSubjectsBetas
                end
            end
        end
    end
end