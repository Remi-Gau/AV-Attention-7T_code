function Vertex_Analysis_07_MVPA_Layer_Plot

clc

StartDirectory = pwd;

SubjectList = [...
    '02';...
%     '03';...
    '04';...
    %     '06';...
    '07';...
    '08';...
%     '09';...
    '11';...
    '12';...
    '13';...
    %     '14';...
    '15';...
    '16'
    ];

Visible = 'off';

FontSize = 6;

% Options for the SVM
opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions


Analysis = struct('name', 'A Stim VS V Stim');
Analysis(end+1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');


ROIs = {...
    'A1'; ... 1
    'TE1.0'; ...
    'TE1.1'; ...
    'TE1.2'; ...
    'PT'; ...
    
    'V1'; ... 6
    'V2'; ...
    'V3d'; ... 8
    'V3v'; ... 9
    'V4d'; ... 10
    'V4v'; ... 11
    'V5';...
    'V3';...
    'V4'};

Layers = [5 7];

%
SaveSufix = '_results_crosslayer';

if opt.fs.do
    SaveSufix = [SaveSufix '_FS'];
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

SaveSufix = [SaveSufix '.mat'];


Sides = {'LeftHS' ; 'RightHS' ; 'BothHS'};


XLabel = char(ROIs);
XLabel = [XLabel(1:5,:);repmat(' ',1,size(XLabel,2));XLabel(6:end,:)];
XLabel = cellstr(XLabel);


% Color for Subjects
COLOR_Subject= [
    0,0,0;
    31,120,180;
    178,223,138;
    %     51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    %     0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;


%%

if ~exist(fullfile(StartDirectory, 'Vertices_Analysis', 'Figures' , 'Data_Block_MVPA_CrossLayers.mat'), 'file')
    
    for SubjInd = 1:size(SubjectList,1)
        
        SubjID = SubjectList(SubjInd,:);
        
        fprintf('\n\nAnalysing subject %s\n', SubjID)
        
        SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
        
        SaveDir = fullfile(SubjectFolder, 'Analysis', 'Vertex', 'SVM');
        
        for iSVM=1:numel(Analysis)
            
            fprintf(' SVM: %s.\n', Analysis(iSVM).name)
            
            for iROI=1:size(ROIs,1)
                
                fprintf('  ROI: %s\n', ROIs{iROI,1})
                
                try
                    
                    load(fullfile(SaveDir, ...
                        ['SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} SaveSufix]), ...
                        'Class_Acc', 'Results');
                    
                    AccROI(:,:,SubjInd,iROI,iSVM) = ...
                        Class_Acc.TotAcc;
                    
                    NbVertices(iROI,SubjInd,1:2) = Results.size;
                    NbVertices(iROI,SubjInd,3) = sum(NbVertices(iROI,SubjInd,1:2));
                    
                    for iLayer = 1:numel(Layers)
                        for iSide = 1:3
                            AccLayer{iLayer}(:,:,SubjInd,iSide,iROI,iSVM) = ...
                                Class_Acc.TotAccLayers{iSide,iLayer};
                            
                        end
                    end
                    clear iLayer iSide
                    
                catch
                    warning('Could not process %s', ...
                        ['SVM_' ROIs{iROI} ': ' Analysis(iSVM).name '_ROI_' ROIs{iROI}])
                    
                    AccROI(:,:,SubjInd,iROI,iSVM) = ...
                        nan(3,2);
                    
                    for iLayer = 1:length(Layers)
                        for iSide = 1:3
                            AccLayer{iLayer}(:,:,SubjInd,iSide,iROI,iSVM) = ...
                                nan(Layers(iLayer));
                        end
                    end
                end
                clear Class_Acc iLayers Results
                
            end
            clear iROI
            
        end
        clear iSVM
        
    end
    clear NbRuns SubjInd
    
    save(fullfile(StartDirectory, 'Vertices_Analysis', 'Figures' , 'Data_Block_MVPA_CrossLayers.mat'))
    
else
    load(fullfile(StartDirectory, 'Vertices_Analysis', 'Figures' , 'Data_Block_MVPA_CrossLayers.mat'))
end

%% Remove subject's values for the ROIs that have less than 80% coverage
Threshold=.8;

SubjectList_Ori = SubjectList;

load(fullfile(StartDirectory,  'ROI_Analysis', 'Figures', '10_layers', 'Data_10_Layers.mat'), 'AllSubjects_Data', 'SubjectList')

Include = [];

for i=[1 3 2 4]
    if i<3
        J=[1 3 4 5 2];
    else
        J=1:7;
    end
    for j=J
        Include(end+1,:) = (AllSubjects_Data(i).SubROI(j).ROI_Coverage(:,1)...
            ./AllSubjects_Data(i).SubROI(j).ROI_Coverage(:,2));
    end
end

A = str2num(SubjectList_Ori);
B = str2num(SubjectList);
LIB = ismember(B,A);
B  = B(LIB,:);
Include = Include(:,LIB);

LIA = ismember(A,B);
temp(:,LIA) = Include;
temp(:,~LIA) = 1;

Include = temp(1:size(temp,1)/2,:);
Include(:,:,2) = temp(1+size(temp,1)/2:end,:);
Include(:,:,3) = mean(Include,3);

clear temp

Include(end+1,:,:) = Include(8,:,:)+Include(9,:,:);
Include(end+1,:,:) = Include(10,:,:)+Include(11,:,:);
Include = Include>Threshold;

SubjectList = SubjectList_Ori;

temp = cat(4,Include,NbVertices>100);
temp = all(temp,4);

Include = temp;

clear A B i j J LOCB Threshold SubjectList_Ori AllSubjects_Data temp


%% Print accuracies
close all

% Color map
X = 0:0.001:1;

R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
clear X


for iSVM=1:numel(Analysis)
    
    for iLayers = 1:length(Layers)
        
        for iSide = 1:size(Sides,1)
            
            temp = Include(:,:,iSide);
            
            % Get data for all SubROIs for that SVM
            AccLayerTmp = squeeze(AccLayer{iLayers}(:,:,:,iSide,:,iSVM));
            AccLayerTmp(AccLayerTmp==0) = NaN;
            
            AccROITmp =  squeeze(AccROI(iSide,iLayers,:,:,iSVM));
            AccROITmp(AccROITmp==0) = NaN;
            
            AllAtChance = squeeze(all(all(AccLayerTmp==0.5,2)))';
            temp(:,:,2) = AllAtChance;
            temp = any(temp,3);
            
            if any(AllAtChance(:))
                warning('An analysis has all its accuracy at chance level.')
                [I,J] = find(AllAtChance);
                for i=1:length(I)
                    fprintf([' ' Sides{iSide} ': subject ' SubjectList(J(i),:) ' for ROI '  ROIs{I(i)} '\n'])
                end
            end
            
            for iROI = 1:size(temp,1)
                if any(~temp(iROI,:))
                    AccLayerTmp(:,:,~temp(iROI,:),iROI) = nan([Layers(iLayers) Layers(iLayers) sum(~temp(iROI,:))]);
                    AccROITmp(~temp(iROI,:),iROI) = NaN;
                end
            end
            

            
            % Compute means accross subjects
            MEAN = squeeze(nanmean(AccLayerTmp,3));
            STD = squeeze(nanstd(AccLayerTmp,3)); %#ok<NASGU>
            SEM = squeeze(nansem(AccLayerTmp,3)); 
            
%             CLIM = [0 1];
%             CLIM2 = [0 1];
            
            CLIM = [min(AccLayerTmp(:)) max(AccLayerTmp(:))];
            CLIM = max(abs(CLIM - .5));
            CLIM = [0.5-CLIM 0.5+CLIM]*10;
            CLIM = [floor(CLIM(1)) ceil(CLIM(2))]/10;
            
            CLIM2 = [min(MEAN(:)) max(MEAN(:))];
            CLIM2 = max(abs(CLIM2 - .5));
            CLIM2 = [0.5-CLIM2 0.5+CLIM2]*10;
            CLIM2 = [floor(CLIM2(1)) ceil(CLIM2(2))]/10;
            
            
            %% Whole ROI accuracies
            figure('name', ['WholeROI_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name],...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
            
            hold on
            grid on
            
            plot([0 13], [0.5 0.5], '--k', 'linewidth', 2)
            
            for SubjInd = 1:size(AccROITmp,1)
                
                plot(1:5,AccROITmp(SubjInd,1:5), 'b', 'linewidth', 1, 'Color', COLOR_Subject(SubjInd,:))
                plot(7:13,AccROITmp(SubjInd,6:12), 'b', 'linewidth', 1, 'Color', COLOR_Subject(SubjInd,:))
                
            end
            
            shadedErrorBar(1:5, nanmean(AccROITmp(:,1:5),1),nansem(AccROITmp(:,1:5),1), ...
                {'LineWidth', 2, 'Color', 'b'}, 1)
            
            shadedErrorBar(7:13, nanmean(AccROITmp(:,6:12),1),nansem(AccROITmp(:,6:12),1), ...
                {'LineWidth', 2, 'Color', 'b'}, 1)
            
            axis([0.5 size(AccROITmp,2)+.5 min(AccROI(:))-.05 1])
            
            set(gca,'tickdir', 'out', ...
                'xtick', 1:size(AccROITmp,2)+1, ...
                'xticklabel', XLabel, ...
                'ytick', 0:.1:1, ...
                'yticklabel', 0:.1:1, 'FontSize', FontSize)
            
            t=ylabel(sprintf([Sides{iSide} '\naccuracy']));
            set(t,'fontsize',FontSize);
            
            print(gcf, fullfile(StartDirectory,  'Vertices_Analysis', 'Figures', ...
                ['WholeROI_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name]), '-dtiff')
            
            
            % Within Layer Accuracies
            figure('name', ['WithinLayer_Audio_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name],...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
            
            iSubPlotWithin=1;
            
            
            WithinLayerPlot(AccLayer, AccLayerTmp, MEAN, SEM, iSubPlotWithin, [1 5], ...
                Layers, iLayers, COLOR_Subject, ROIs, FontSize);
            
            subplot(1, length(1:5), 1)
            t=ylabel(sprintf([Sides{iSide} '\nWithin ROI accuracy']));
            set(t,'fontsize',FontSize);
            
            
            print(gcf, fullfile(StartDirectory,  'Vertices_Analysis', 'Figures', ...
                ['WithinLayer_Audio_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name '.tif']), '-dtiff')
            
            
            
            figure('name', ['WithinLayer_Visual_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name],...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
            
            iSubPlotWithin=1;
            
            WithinLayerPlot(AccLayer, AccLayerTmp, MEAN, SEM, iSubPlotWithin, [6 12], ...
                Layers, iLayers, COLOR_Subject, ROIs, FontSize);
            
            
            subplot(1, length(6:12), 1)
            t=ylabel(sprintf([Sides{iSide} '\nWithin ROI accuracy']));
            set(t,'fontsize',FontSize);
            
            
            print(gcf, fullfile(StartDirectory,  'Vertices_Analysis', 'Figures', ...
                ['WithinLayer_Visual_' Sides{iSide} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name '.tif']), '-dtiff')
            
            
            
            %% Cross Layer Accuracies
            figure('name', ['CrossLayer_' Sides{iSide,1} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name],...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible)
            
            iSubPlot=1;
            
            % Plot each subjects
            for SubjInd = 1:size(AccLayer{1},3)
                for iROI = 1:size(ROIs,1)
                    
                    subplot(2+size(AccLayer{1},3), 1+size(AccLayer{1},5), iSubPlot)
                    
                    colormap([R' G' B']);
                    
                    imagesc(flipud(AccLayerTmp(:,:,SubjInd,iROI)), CLIM)
                    
                    axis square
                    
                    set(gca,'tickdir', 'out', ...
                        'xtick', 1:Layers(iLayers), ...
                        'xticklabel', [], ...
                        'ytick', 1:Layers(iLayers), ...
                        'yticklabel', [], 'FontSize', FontSize)
                    
                    if iROI==1
                        t=ylabel(sprintf('Subj %s', SubjectList(SubjInd,:)));
                        set(t,'fontsize',FontSize+2);
                    end
                    
                    iSubPlot=iSubPlot+1;
                    
                    if iROI==5
                        if SubjInd==size(AccLayer{1},3)
                        subplot(2+size(AccLayer{1},3), 1+size(AccLayer{1},5), iSubPlot)
                        
                        
                        colormap([R' G' B']);
                        imagesc(repmat([1000:-1:1]',1,100)) %#ok<NBRAK>
                        axis square
                        
                        set(gca,'tickdir', 'out', ...
                            'ytick', [1 500 1000], ...
                            'yticklabel', [CLIM(2) 0.5 CLIM(1)], ...
                            'xtick', [], ...
                            'xticklabel', [], 'FontSize', FontSize)
                        end
                        
                        iSubPlot=iSubPlot+1;
                    end
                end
            end
            
            % Skip a line
            iSubPlot=iSubPlot+1+size(AccLayer{1},5);
            
            % Plot means
            for iROI = 1:size(ROIs,1)
                
                subplot(2+size(AccLayer{1},3), 1+size(AccLayer{1},5), iSubPlot)
                
                colormap([R' G' B']);
                
                imagesc(flipud(MEAN(:,:,iROI)), CLIM2)
                
                axis square
                
                t=xlabel(sprintf('Test layer\n\n %s', ROIs{iROI}));
                set(t,'fontsize',FontSize+2);
                
                if iROI==1
                    t=ylabel(sprintf('MEAN\n\nTraining layer'));
                    set(t,'fontsize',FontSize+2);
                end
                
                set(gca,'tickdir', 'out', ...
                    'xtick', 1:2:Layers(iLayers), ...
                    'xticklabel', 1:2:Layers(iLayers), ...
                    'ytick', 1:2:Layers(iLayers), ...
                    'yticklabel', Layers(iLayers):-2:1, 'FontSize', FontSize)
                
                iSubPlot=iSubPlot+1;
                
                if iROI==5
                    subplot(2+size(AccLayer{1},3), 1+size(AccLayer{1},5), iSubPlot)
                    
                    
                    colormap([R' G' B']);
                    imagesc(repmat([1000:-1:1]',1,100)) %#ok<NBRAK>
                    axis square
                    
                    set(gca,'tickdir', 'out', ...
                        'ytick', [1 500 1000], ...
                        'yticklabel', [CLIM2(2) 0.5 CLIM2(1)], ...
                        'xtick', [], ...
                        'xticklabel', [], 'FontSize', FontSize)
                    
                    iSubPlot=iSubPlot+1;
                end
            end
            
            
            print(gcf, fullfile(StartDirectory,  'Vertices_Analysis', 'Figures', ['CrossLayer_' Sides{iSide,1} '_' num2str(Layers(iLayers)) '_Layers_' Analysis(iSVM).name '.tif']), '-dtiff')
            
            
        end
    end
end


%% Color map
figure('name', 'Color Map', 'Position', [400, 400, 100, 400], ...
    'Color', [1 1 1], 'Visible', 'on')

colormap([R' G' B']);

imagesc(repmat([1000:-1:1]',1,100)) %#ok<NBRAK>

set(gca,'tickdir', 'out', ...
    'ytick', 1:100:1000, ...
    'yticklabel', 1:-0.1:0, ...
    'xtick', [], ...
    'xticklabel', [])

print(gcf, fullfile(pwd,  'Vertices_Analysis', 'Figures', 'Legend.tif'), '-dtiff')

end

function [iSubPlot] = WithinLayerPlot(AccLayer, AllAcc, Mean, Error, iSubPlot, ROI_Nb, Layers, iLayers, COLOR_Subject, ROIs, FontSize)

ROI_Nb = ROI_Nb(1):ROI_Nb(2);

for iROI = 1:length(ROI_Nb)
    
    subplot(1, length(ROI_Nb), iSubPlot)
    
    hold on
    
    shadedErrorBar(1:Layers(iLayers), diag(Mean(:,:,ROI_Nb(iROI))),...
        flipud(diag(Error(:,:,ROI_Nb(iROI)))), ...
        {'LineWidth', 2, 'Color', 'b'}, 0)
    
    plot([1 Layers(iLayers)], [0.5 0.5], '--k', 'linewidth', 1)
    
    for SubjInd = 1:size(AccLayer{1,1,1},3)
        
        plot(diag(AllAcc(:,:,SubjInd,ROI_Nb(iROI))), 'Color', COLOR_Subject(SubjInd,:), ...
            'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
            'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 2)
    end
    
    
    axis([0.8 Layers(iLayers)+.2  min(AllAcc(:))-.05 1])
    
    if iSubPlot==1
        set(gca,'tickdir', 'out', ...
            'xtick', 1:Layers(iLayers), ...
            'xticklabel', 1:Layers(iLayers), ...
            'ytick', 0:.1:1, ...
            'yticklabel', 0:.1:1, 'FontSize', FontSize)
    else
        set(gca,'tickdir', 'out', ...
            'xtick', 1:Layers(iLayers), ...
            'xticklabel', [], ...
            'ytick', 0:.1:1, ...
            'yticklabel',[], 'FontSize', FontSize)
    end
    
    t=xlabel(sprintf('%s', ROIs{ROI_Nb(iROI)}));
    set(t,'fontsize',FontSize+2);
    
    iSubPlot=iSubPlot+1;
    
end
end