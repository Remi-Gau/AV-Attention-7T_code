clear
clc

StartDirectory = fullfile(pwd, '..', '..', '..');

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

for iSubj=1:size(SubjectList,1)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
% [LIA,LOCB] = ismember( ones(1,size(SubjectList,1)),ToPermute, 'rows');
% ToPermute(LOCB,:) = [];


FigDim = [100, 100, 1300, 800];

FontSize = 6;

% Options for the SVM
opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

opt.scaling.idpdt = 1;
opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;

% ROIs = {...
%     'V1'; ...
%     'V2-3'; ...
%     'A1'; ...
%     'PT'; ...
%     };

ROIs = {...
    'A1'; ...
    'PT'; ...
    'V1'; ...
    'V2-3'; ...
    'V1_act'; ...
    'V1_deact'; ...
    'V23_act'; ...
    'V23_deact'; ...
    };

Analysis(1) = struct('name', 'A Stim VS AV Stim', 'ROIs', [1 2]); %#ok<*AGROW>
Analysis(end+1) = struct('name', 'V Stim VS AV Stim', 'ROIs', 3:8);
Analysis(end+1) = struct('name', 'A Att VS V Att', 'ROIs', 1:8);
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)', 'ROIs', 1:8);
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)', 'ROIs', 1:8);
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)', 'ROIs', 1:8);

NbLayers = 7;
NbLayers = NbLayers+1;

LayersDepth = round(linspace(0,100,NbLayers));
LayersDepth([1 end]) = [];

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

SaveSufix = [SaveSufix '_FWHM_' FFX{1} '_Layers_' num2str(NbLayers) '.mat'];


%%

if ~exist(fullfile(StartDirectory, 'Figures' , 'Data_MVPA_AVT_CrossLayers.mat'), 'file')
    
    for SubjInd = 1:size(SubjectList,1)
        
        SubjID = SubjectList(SubjInd,:);
        
        fprintf('\n\nAnalysing subject %s\n', SubjID)
        
        SubjectFolder = fullfile(StartDirectory,...
            'Subjects_Data', ['Subject_' SubjID]);
        
        SaveDir = fullfile(SubjectFolder,'Transfer', 'SVM');
        
        for iSVM=1:numel(Analysis)
            
            fprintf(' SVM: %s.\n', Analysis(iSVM).name)
            
            for iROI=1:size(ROIs,1)
                
                fprintf('  ROI: %s\n', ROIs{iROI,1})
                
                clear CV Results temp
                
                try
                    
                    load(fullfile(SaveDir, ...
                        ['SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} SaveSufix]), 'Results');
                    
                    CV = Results.session(end).rand.perm.CV;
                    Class_Acc(:,iSVM,iROI,SubjInd) = nanmean([CV(:).acc]);
                    for iCV=1:numel(CV)
                        temp(:,:,iCV) = CV(iCV).layers.acc(2:end-1,2:end-1);
                    end
                    
                    if  all(all(nanmean(temp(:,:,:),3)==0.5*ones(NbLayers-2)))
                        Class_Acc(:,iSVM,iROI,SubjInd) = nan(5,1);
                        Class_Acc_Layers(:,:,iSVM,iROI,SubjInd) = nan(NbLayers-2);
                    else
                        Class_Acc_Layers(:,:,iSVM,iROI,SubjInd) = nanmean(temp,3);
                        Class_Acc(2:7,iSVM,iROI,SubjInd) = diag(Class_Acc_Layers(:,:,iSVM,iROI,SubjInd));
                    end
                    
                    NbVertices(iROI,SubjInd,1:2) = Results.size;
                    NbVertices(iROI,SubjInd,3) = sum(NbVertices(iROI,SubjInd,1:2));
                    
                    
                catch
                    warning('Could not process %s', ...
                        ['SVM_' ROIs{iROI} ': ' Analysis(iSVM).name '_ROI_' ROIs{iROI}])
                    
                    Class_Acc(:,iSVM,iROI,SubjInd) = nan(5,1);
                    Class_Acc_Layers(:,:,iSVM,iROI,SubjInd) = nan(4);
                    
                end
                clear iLayers
                
            end
            clear iROI
            
        end
        clear iSVM
        
    end
    clear NbRuns SubjInd
    
    save(fullfile(StartDirectory, 'Figures' , 'Data_MVPA_AVT_CrossLayers.mat'))
    
else
    load(fullfile(StartDirectory, 'Figures' , 'Data_MVPA_AVT_CrossLayers.mat'), 'Class_Acc', 'Class_Acc_Layers')
end


%% Print accuracies
close all

% Color map
X = 0:0.001:1;

R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
ColorMap1 = [R' G' B'];

R = 1 - 0.392*(1 + erf((X - 0.869)/ 0.255));
G = 1.021 - 0.456*(1 + erf((X - 0.527)/ 0.376));
B = 1 - 0.493*(1 + erf((X - 0.272)/ 0.309));
ColorMap2 = [R' G' B'];
ColorMap2 = flipud(ColorMap2);
clear X

Visible = 'on';


for iSVM=1:numel(Analysis)
    
    close all
    
    % Get data for all SubROIs for that SVM
    AccLayerTmp = squeeze(Class_Acc_Layers(:,:,iSVM,:,:));
    AccLayerTmp(AccLayerTmp==0) = NaN;
    
    AccROITmp =  squeeze(Class_Acc(:,iSVM,:,:));
    AccROITmp(AccROITmp==0) = NaN;
    
    % Compute means accross subjects
    MEAN = squeeze(nanmean(AccLayerTmp,4));
    STD = squeeze(nanstd(AccLayerTmp,4)); %#ok<NASGU>
    SEM = squeeze(nansem(AccLayerTmp,4));
    
    %% Cross Layer Accuracies
    
    for iROI = Analysis(iSVM).ROIs
        
        
        %% Plot all subjects
        figure('name', ['All_Subjects - CrossLayer : ' Analysis(iSVM).name ' - ' ROIs{iROI}],...
            'Position', FigDim, 'Color', [1 1 1], 'Visible', Visible)
        
        iSubPlot=1;
        
        tmp = abs(AccLayerTmp(:,:,iROI,:));
        MAX = round((max(tmp(:))-.5)*10)/10;
        CLIM = [0.5-MAX 0.5+MAX];
        
        for SubjInd = 1:size(AccLayerTmp,4)
            
            subplot(3, 4, iSubPlot)
            
            colormap(ColorMap1);
            
            imagesc(flipud(AccLayerTmp(:,:,iROI,SubjInd)), CLIM)
            
            axis square
            
            set(gca,'tickdir', 'out', ...
                'xtick', 1:NbLayers-2, ...
                'xticklabel', LayersDepth, ... 
                'ytick', 1:NbLayers-2, ...
                'yticklabel', fliplr(LayersDepth), 'FontSize', FontSize)
            
            t=title(sprintf('Subj %s', SubjectList(SubjInd,:)));
            set(t,'fontsize',FontSize+2);
            
            t=xlabel('Test depth');
            set(t,'fontsize',FontSize+2);
            
            t=ylabel(sprintf('Training depth'));
            set(t,'fontsize',FontSize+2);
            
            iSubPlot=iSubPlot+1;
            
        end
        
        subplot(3, 4, iSubPlot)
        
        colormap(ColorMap1);
        imagesc(repmat([CLIM(2):-.01:CLIM(1)]', [1,200]), CLIM)
        
        axis square
        
        set(gca,'tickdir', 'out', 'xtick', [],'xticklabel',  [], ...
            'ytick', linspace(1,numel(CLIM(2):-.01:CLIM(1)),5),...
            'yticklabel', linspace(CLIM(2),CLIM(1),5), ...
            'ticklength', [0.01 0.01], 'fontsize', 8, 'YAxisLocation','right')
        
        box off
        
        
        mtit([Analysis(iSVM).name ' - ' strrep(ROIs{iROI},'_', ' ')], 'xoff', 0, 'yoff', +0.03)
        
        Name = fullfile(StartDirectory, 'Figures',...
            ['AllSubjects_' strrep(Analysis(iSVM).name, ' ', '_') '_' ROIs{iROI} '_Surfaces_CrossLayers.tif']);
        Name = strrep(Name, '(', '-');
        Name = strrep(Name, ')', '');
        print(gcf, Name, '-dtiff')
        
        
        %% Plot MEAN
        figure('name', ['MEAN - CrossLayer : ' Analysis(iSVM).name ' - ' ROIs{iROI}],...
            'Position', FigDim, 'Color', [1 1 1], 'Visible', Visible)
        
        tmp = abs(MEAN(:,:,Analysis(iSVM).ROIs));
%         CLIM = [floor(100*min(tmp(:)))/100 floor(100*max(tmp(:)))/100];
        CLIM = [0.5 floor(100*max(tmp(:)))/100];
        
        
        
        colormap(ColorMap2);
        imagesc(flipud(MEAN(:,:,iROI)), CLIM)
        
        tmp = squeeze(AccLayerTmp(:,:,iROI,:));
        Pos = (NbLayers-2):-1:1;

        for iPerm = 1:size(ToPermute,1)
%             tmp2(1,1,:) = ToPermute(iPerm,:);
%             Perms(:,:,iPerm) = mean((tmp-.5).*repmat(tmp2,[NbLayers-2,NbLayers-2,1]),3); 
            tmp1(1,1,:) = ToPermute(iPerm,:);
            tmp2 = mean((tmp-.5).*repmat(tmp1,[NbLayers-2,NbLayers-2,1]),3);
            Perms(iPerm,1) = max(abs(tmp2(:))); %#ok<*SAGROW>
        end
        
        for x=1:(NbLayers-2)
            for y=1:(NbLayers-2)
                
%               [H,P] = ttest(squeeze(tmp(x,y,:)), .5, 'alpha', .05/(NbLayers-2)^2);
%               [H,P] = ttest(squeeze(tmp(x,y,:)), .5, 'alpha', .05);
              
%               P = sum(Perms(x,y,:)>(mean(tmp(x,y,:))-.5))/numel(Perms(x,y,:));

                P = sum(Perms>(mean(tmp(x,y,:))-.5))/numel(Perms);

                if P<0.001
                    Sig = sprintf('p<0.001');
                else
                    Sig = sprintf('p=%.3f',P);
                end
                t = text(x-.15,Pos(y),sprintf(Sig));
                set(t,'fontsize',8);
                
                if P<(.05)
                    set(t,'color','g');
                end

            end
        end
        
        axis square
        
        t=xlabel('Test depth');
        set(t,'fontsize',12);
        
        t=ylabel(sprintf('Training depth'));
        set(t,'fontsize',12);
        
        set(gca,'tickdir', 'out', ...
            'xtick', 1:NbLayers-2, ...
            'xticklabel', LayersDepth, ...
            'ytick', 1:NbLayers-2, ...
            'yticklabel', fliplr(LayersDepth), 'FontSize', 12)
        
        
        pause(.5)
        
        ax=gca;
        axPos = ax.Position;
        axPos(1) = axPos(1)+axPos(3)-0.05;
        axPos(3) = .05;
        axes('Position',axPos);
        
        colormap(ColorMap2);
        imagesc(repmat([CLIM(2):-.001:CLIM(1)]', [1,200]), CLIM)
        set(gca,'tickdir', 'out', 'xtick', [],'xticklabel',  [], ...
            'ytick', linspace(1,numel(CLIM(2):-.001:CLIM(1)),5),...
            'yticklabel', linspace(CLIM(2),CLIM(1),5), ...
            'ticklength', [0.01 0.01], 'fontsize', 8, 'YAxisLocation','right')
        box off
        
        
        mtit(['Mean - ' Analysis(iSVM).name ' - ' strrep(ROIs{iROI},'_', ' ')], 'xoff', 0, 'yoff', +0.03)
        
        Name = fullfile(StartDirectory, 'Figures',...
            ['Mean_' strrep(Analysis(iSVM).name, ' ', '_') '_' ROIs{iROI} '_Surfaces_CrossLayers.tif']);
        Name = strrep(Name, '(', '-');
        Name = strrep(Name, ')', '');
        print(gcf, Name, '-dtiff')
        
        
        
    end
    
    clear AccLayerTmp MEAN AccROITmp
end