function Analysis_MVPA_03_MVPA_Vertices_Plot

clc

StartDirectory = fullfile(pwd, '..', '..', '..', '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %'14';...
    '15';...
    '16'
    ];

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



Analysis(1) = struct('name', 'A Stim VS AV Stim'); %#ok<*AGROW>
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');



ROIs = {...
    'V1'; ...
    'V2-3'; ...
    'A1'; ...
    'PT'; ...
    };

NbLayers = 7;
NbLayers = NbLayers+1;

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

% %% Remove subject's values for the ROIs that have less than 80% coverage
% Threshold=.8;
%
% SubjectList_Ori = SubjectList;
%
% load(fullfile(StartDirectory,  'ROI_Analysis', 'Figures', '10_layers', 'Data_10_Layers.mat'), 'AllSubjects_Data', 'SubjectList')
%
% Include = [];
%
% for i=[1 3 2 4]
%     if i<3
%         J=[1 3 4 5 2];
%     else
%         J=1:7;
%     end
%     for j=J
%         Include(end+1,:) = (AllSubjects_Data(i).SubROI(j).ROI_Coverage(:,1)...
%             ./AllSubjects_Data(i).SubROI(j).ROI_Coverage(:,2));
%     end
% end
%
% A = str2num(SubjectList_Ori);
% B = str2num(SubjectList);
% LIB = ismember(B,A);
% B  = B(LIB,:);
% Include = Include(:,LIB);
%
% LIA = ismember(A,B);
% temp(:,LIA) = Include;
% temp(:,~LIA) = 1;
%
% Include = temp(1:size(temp,1)/2,:);
% Include(:,:,2) = temp(1+size(temp,1)/2:end,:);
% Include(:,:,3) = mean(Include,3);
%
% clear temp
%
% Include(end+1,:,:) = Include(8,:,:)+Include(9,:,:);
% Include(end+1,:,:) = Include(10,:,:)+Include(11,:,:);
% Include = Include>Threshold;
%
% SubjectList = SubjectList_Ori;
%
% temp = cat(4,Include,NbVertices>100);
% temp = all(temp,4);
%
% Include = temp;
%
% clear A B i j J LOCB Threshold SubjectList_Ori AllSubjects_Data temp


%% Print accuracies
close all

% Color map
X = 0:0.001:1;

R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
clear X

% Color for Subjects
COLOR_Subject= [
    0,0,0;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

Visible = 'off';

Transparent = 1;

for iSVM=1:numel(Analysis)
    
    
    %             temp = Include(:,:,iSide);
    
    % Get data for all SubROIs for that SVM
    AccLayerTmp = squeeze(Class_Acc_Layers(:,:,iSVM,:,:));
    AccLayerTmp(AccLayerTmp==0) = NaN;
    
    AccROITmp =  squeeze(Class_Acc(:,iSVM,:,:));
    AccROITmp(AccROITmp==0) = NaN;
    
    % Compute means accross subjects
    MEAN = squeeze(nanmean(AccLayerTmp,4));
    STD = squeeze(nanstd(AccLayerTmp,4)); %#ok<NASGU>
    SEM = squeeze(nansem(AccLayerTmp,4));
    
    
    
    CLIM = [0 1];
    CLIM2 = [0 1];
    
%     CLIM = [min(AccLayerTmp(:)) max(AccLayerTmp(:))];
%     CLIM = max(abs(CLIM - .5));
%     CLIM = [0.5-CLIM 0.5+CLIM]*10;
%     CLIM = [floor(CLIM(1)) ceil(CLIM(2))]/10;
%     
%     CLIM2 = [min(MEAN(:)) max(MEAN(:))];
%     CLIM2 = max(abs(CLIM2 - .5));
%     CLIM2 = [0.5-CLIM2 0.5+CLIM2]*10;
%     CLIM2 = [floor(CLIM2(1)) ceil(CLIM2(2))]/10;
    
    
    %%
    figure('position', FigDim, 'name', ['All_Subjects : ' Analysis(iSVM).name], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    mn = length(ROIs);
    n  = round(mn^0.4);
    m  = ceil(mn/n);
    
    iSubPlot=1;
    
    for FWHM = 1
        
        for iROI = 1:numel(ROIs)
            
            subplot(n,m,iSubPlot)
            
            hold on
            grid on
            
            for SubjInd = 1:size(AccROITmp,3)
                plot(1:6, (AccROITmp(2:end,iROI,SubjInd)), ...
                    '--','marker', '.', 'color', COLOR_Subject(SubjInd,:), ...
                    'markerfacecolor', COLOR_Subject(SubjInd,:),'linewidth', .5)
            end
            
            plot([1 6], [0.5 0.5], '--k', 'linewidth', 2)
            
            axis([0.5 6.5 0 1])
            
            set(gca,'tickdir', 'out', ...
                'xtick', 1:6, ...
                'xticklabel', 1:6, ...
                'ytick', 0:.1:1, ...
                'yticklabel', 0:.1:1,...
                'fontsize', 8)
            
            title(ROIs{iROI},'fontsize', 10)
            
            shadedErrorBar(1:6, mean(squeeze(AccROITmp(2:end,iROI,:)),2),...
                nansem(squeeze(AccROITmp(2:end,iROI,:)),2), ...
                {'LineWidth', 3, 'Color', [0 0 1]}, Transparent)
            
            iSubPlot = iSubPlot + 1;
        end
        
    end
    
    mtit(Analysis(iSVM).name, 'xoff', 0, 'yoff', +0.02)
    
    Name = fullfile(StartDirectory, 'Figures',...
        ['AllSubjects_' strrep(Analysis(iSVM).name, ' ', '_') '_Surfaces.pdf']);
    Name = strrep(Name, '(', '-');
    Name = strrep(Name, ')', '');
    print(gcf, Name, '-dpdf')
    
    
    
    %% Cross Layer Accuracies
    figure('name', ['All_Subjects - CrossLayer : ' Analysis(iSVM).name],...
        'Position', FigDim, 'Color', [1 1 1], 'Visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    iSubPlot=1;
    
    % Plot each subjects
    for SubjInd = 1:size(AccLayerTmp,4)
        for iROI = 1:size(ROIs,1)
            
            subplot(2+size(AccLayerTmp,4), size(ROIs,1), iSubPlot)
            
            colormap([R' G' B']);
            
            imagesc(flipud(AccLayerTmp(:,:,iROI,SubjInd)), CLIM)
            
            axis square
            
            set(gca,'tickdir', 'out', ...
                'xtick', [], ...
                'xticklabel', [], ...
                'ytick', [], ...
                'yticklabel', [], 'FontSize', FontSize-4)
            
            if iROI==1
                t=ylabel(sprintf('Subj %s', SubjectList(SubjInd,:)));
                set(t,'fontsize',FontSize+2);
            end
            
            iSubPlot=iSubPlot+1;
            
        end
        
    end
    
    % Skip a line
    iSubPlot=iSubPlot+size(ROIs,1);
    
    % Plot means
    for iROI = 1:size(ROIs,1)
        
        subplot(2+size(AccLayerTmp,4), size(ROIs,1), iSubPlot)
        
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
            'xtick', 1:NbLayers-2, ...
            'xticklabel', 1:NbLayers-2, ...
            'ytick', 1:NbLayers-2, ...
            'yticklabel', NbLayers-2:-1:1, 'FontSize', FontSize)
        
        iSubPlot=iSubPlot+1;
        
    end
    
    mtit(Analysis(iSVM).name, 'xoff', 0, 'yoff', +0.02)
    
    Name = fullfile(StartDirectory, 'Figures',...
        ['AllSubjects_' strrep(Analysis(iSVM).name, ' ', '_') '_Surfaces_CrossLayers.pdf']);
    Name = strrep(Name, '(', '-');
    Name = strrep(Name, ')', '');
    print(gcf, Name, '-dpdf')
    
    
    clear AccLayerTmp MEAN AccROITmp
    
end

Command =[];

for iSVM=1:numel(Analysis)
    
    Name = fullfile(StartDirectory, ...
        [strrep(Analysis(iSVM).name, ' ', '_') '.pdf']);
    Name = strrep(Name, '(', '-');
    Name = strrep(Name, ')', '');
    Command = [Command ' ' Name];
    
    Name = fullfile(StartDirectory, 'Figures',...
        ['AllSubjects_' strrep(Analysis(iSVM).name, ' ', '_') '_Surfaces.pdf']);
    Name = strrep(Name, '(', '-');
    Name = strrep(Name, ')', '');
    Command = [Command ' ' Name];
    
    Name = fullfile(StartDirectory, 'Figures',...
        ['AllSubjects_' strrep(Analysis(iSVM).name, ' ', '_') '_Surfaces_CrossLayers.pdf']);
    Name = strrep(Name, '(', '-');
    Name = strrep(Name, ')', '');
    Command = [Command ' ' Name];
    
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, 'AllSubjects_AV_MVPA.pdf') ' ' Command])

% %% Color map
figure('name', 'Color Map', 'Position', [400, 400, 100, 400], ...
    'Color', [1 1 1], 'Visible', 'on')

colormap([R' G' B']);

imagesc(repmat([1000:-1:1]',1,100)) %#ok<NBRAK>

set(gca,'tickdir', 'out', ...
    'ytick', 1:100:1000, ...
    'yticklabel', 1:-0.1:0, ...
    'xtick', [], ...
    'xticklabel', [])

% print(gcf, fullfile(pwd, 'Legend.tif'), '-dtiff')

end