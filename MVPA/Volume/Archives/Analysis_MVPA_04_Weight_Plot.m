%%
clc; clear; close all;

StartDirectory = pwd;

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %     '06';...
    %     '07';... % Re run: prob with V4d L
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %     '14';...
    '15';...
    '16'
    ];


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


Analysis = struct('name', 'A Stim VS V Stim');
Analysis(end+1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');


ROIs = { ...
    'A_lh', 1:5; ...
    'V_lh', 6:12; ...
    'A_rh', 13:17; ...
    'V_rh', 18:24; ...
    };


SubROIs = {...
    'A1L'; ... 1
    'TE1.0L'; ...
    'TE1.1L'; ...
    'TE1.2L'; ...
    'PTL'; ...
    
    'V1L'; ... 6
    'V2L'; ...
    'V3dL'; ...
    'V3vL'; ...
    'V4dL'; ...
    'V4vL'; ...
    'V5L'; ...
    
    'A1R'; ... 13
    'TE1.0R'; ...
    'TE1.1R'; ...
    'TE1.2R'; ...
    'PTR'; ...
    
    'V1R'; ... 18
    'V2R'; ...
    'V3dR'; ...
    'V3vR'; ...
    'V4dR'; ...
    'V4vR'; ...
    'V5vR'};


Data(1) = struct('MEAN',[],'STD',[],'SEM',[]);
Data(2) = struct('MEAN',[],'STD',[],'SEM',[]);
Data(3) = struct('MEAN',[],'STD',[],'SEM',[]);



%%
for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n\nAnalysing subject %s\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    SaveDir = fullfile(SubjectFolder, 'Analysis', 'ROI');
    
    cd(SaveDir)
    %     delete SVM_Weight.mat
    delete SVM_Restrict_Weight.mat
    
    if exist('SVM_Weight.mat', 'file')
        load SVM_Weight.mat
        
    else
        
        WeightMean{1,1} = nan(5, size(SubROIs,1),numel(Analysis), 2);
        WeightStd{1,1} = nan(5, size(SubROIs,1),numel(Analysis), 2);
        WeightSEM{1,1} = nan(5, size(SubROIs,1),numel(Analysis), 2);
        WeightMean{2,1} = nan(7, size(SubROIs,1),numel(Analysis), 2);
        WeightStd{2,1} = nan(7, size(SubROIs,1),numel(Analysis), 2);
        WeightSEM{2,1} = nan(7, size(SubROIs,1),numel(Analysis), 2);
        WeightMean{3,1} = nan(11, size(SubROIs,1),numel(Analysis), 2);
        WeightStd{3,1} = nan(11, size(SubROIs,1),numel(Analysis), 2);
        WeightSEM{3,1} = nan(11, size(SubROIs,1),numel(Analysis), 2);
        
        for iSVM=1:numel(Analysis)
            %%
            
            close all
            
            fprintf(' SVM: %s.\n', Analysis(iSVM).name)
            
            for iROI=1:size(ROIs,1)
                
                for iSubROI = 1:length(ROIs{iROI,2})
                    
                    fprintf('  ROI: %s\n', SubROIs{ROIs{iROI,2}(iSubROI),1})
                    
                    Layers(:,:,:,1) = spm_read_vols(spm_vol(fullfile(SaveDir, ROIs{iROI}, 'T1_04_Layers.nii')));
                    Layers(:,:,:,2) = spm_read_vols(spm_vol(fullfile(SaveDir, ROIs{iROI}, 'T1_06_Layers.nii')));
                    Layers(:,:,:,3) = spm_read_vols(spm_vol(fullfile(SaveDir, ROIs{iROI}, 'T1_10_Layers.nii')));
                    
                    Weights = spm_read_vols(spm_vol(...
                        fullfile(SaveDir, ROIs{iROI}, 'SVM', ...
                        [ROIs{iROI} ': ' Analysis(iSVM).name '_' SubROIs{ROIs{iROI,2}(iSubROI),1} ...
                        '_weight.nii'])));
                    
                    Weights(:,:,:,2) = spm_read_vols(spm_vol(...
                        fullfile(SaveDir, ROIs{iROI}, 'SVM', ...
                        [ROIs{iROI} ': ' Analysis(iSVM).name '_' SubROIs{ROIs{iROI,2}(iSubROI),1} ...
                        '_Restrict_weight.nii'])));
                    
                    iSubplot = 1;
                    
                    figure('Name', [Analysis(iSVM).name ' / ' SubROIs{ROIs{iROI,2}(iSubROI),1} ' : weights distribution'], ...
                        'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible); %#ok<UNRCH>
                    
                    for WeightVol = 1:size(Weights,4)
                        
                        tmp2 = Weights(:,:,:,WeightVol);
                        
                        for iLayerVol = 1:size(Layers,4)
                            
                            LayerVol=Layers(:,:,:,iLayerVol);
                            NbLayers = max(LayerVol(:));
                            
                            Values2Plot={};
                            
                            for iLayer=NbLayers:-1:0
                                
                                tmp = LayerVol==iLayer;
                                tmp(:,:,:,2) = logical(tmp2);
                                tmp = all(tmp,4);
                                
                                WeightMean{iLayerVol,1}(iLayer+1, ROIs{iROI,2}(iSubROI),iSVM, WeightVol) = mean(tmp2(tmp)); %#ok<*NASGU>
                                WeightStd{iLayerVol,1}(iLayer+1, ROIs{iROI,2}(iSubROI),iSVM, WeightVol) = std(tmp2(tmp));
                                WeightSEM{iLayerVol,1}(iLayer+1, ROIs{iROI,2}(iSubROI),iSVM, WeightVol) = nansem(tmp2(tmp));
                                
                                Values2Plot{end+1,1} = tmp2(tmp); %#ok<*UNRCH>
                                
                            end
                            clear iLayer
                            
                            subplot(size(Weights,4),size(Layers,4),iSubplot)
                            
                            grid on
                            hold on
                            
                            distributionPlot(Values2Plot, 'showMM', 1, 'histOpt', 0, 'divFactor', 70, ...
                                'globalNorm', 0, 'color', 'b')
                            
                            set(gca,'tickdir', 'out', 'xtick', 1:NbLayers+1 ,'xticklabel', NbLayers:-1:0, ...
                                'ticklength', [0.01 0.01], 'fontsize', 10)
                            
                            t=ylabel('Weight');
                            set(t,'fontsize', 10);
                            
                            t=xlabel('Layers');
                            set(t,'fontsize', 10);
                            
                            clear Values2Plot
                            
                            iSubplot = iSubplot+1;
                            
                        end
                        clear iLayerVol
                        
                    end
                    clear WeightVol
                    
                    print(gcf, fullfile(SaveDir, ROIs{iROI}, 'SVM', ...
                        [Analysis(iSVM).name '_ROI_' SubROIs{ROIs{iROI,2}(iSubROI),1} '.tif']), '-dtiff')
                    
                end
                clear Weights iSubROI Layers
                
            end
            clear iROI
            
        end
        clear iSVM
        
        save(fullfile(SaveDir, 'SVM_Weight.mat'), 'ROIs', 'SubROIs', 'WeightMean', ...
            'WeightStd', 'WeightSEM')
        
    end
    
    
    for i=1:3
        Data(i).MEAN(:,:,:,:,end+1) = WeightMean{i,1};
        Data(i).STD(:,:,:,:,end+1) = WeightStd{i,1};
        Data(i).SEM(:,:,:,:,end+1) = WeightSEM{i,1};
    end
    
    
    
end
clear iSVM



%%
Visible = 'off';

FontSize = 6;

for iSVM = 1:numel(Analysis)
    
    for iLayers = 1:3
        
        for iRestrict = 1:2
            
            clear DATA2PLOT
            
            DATA2PLOT = squeeze(Data(iLayers).MEAN(:,:,iSVM,iRestrict,:));
            
            if iRestrict==1
                FigName = ['SVM_Weight_' num2str(size(DATA2PLOT,1)-1) '_Layers_' Analysis(iSVM).name];
            else
                FigName = ['SVM_Weight_Restrict_' num2str(size(DATA2PLOT,1)-1) '_Layers_' Analysis(iSVM).name];
            end
            
            figure('name', FigName,...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible)
            
            
            for iSubROI = 1:size(DATA2PLOT,2)
                
                tmp = squeeze(DATA2PLOT(:,iSubROI,:));
                
                tmp(tmp==0)=NaN;
                
                tmp=fliplr(tmp);
                
                %%
                subplot(2,size(DATA2PLOT,2)/2,iSubROI)
                
                hold on
                grid on
                
                errorbar(0.9, nanmean(tmp(:)), nanstd(tmp(:)), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 3)
                
                shadedErrorBar(2:size(tmp,1), nanmean(tmp(1:end-1,:),2),nansem(tmp(1:end-1,:),2), ...
                    {'LineWidth', 3, 'Color', 'b'}, 1)
                
                for SubjInd = 1:size(tmp,2)
                    plot(1.5,nanmean(tmp(:, SubjInd)), ' o', ...
                        'Color', COLOR_Subject(SubjInd,:), 'MarkerFaceColor', COLOR_Subject(SubjInd,:),...
                        'MarkerSize', 2)
                    plot(2:size(tmp,1),tmp(1:end-1, SubjInd), 'b', 'linewidth', 1.5, 'Color', COLOR_Subject(SubjInd,:))
                end
                
                axis([0.5 size(tmp,1) 0 max(tmp(:))])
                
                set(gca,'tickdir', 'out', ...
                    'xtick', [2 3 size(tmp,1)], ...
                    'xticklabel', ['O';'W';'P'],'fontsize',FontSize)
                
                t=xlabel(SubROIs{iSubROI}(1:end-1));
                set(t,'fontsize',FontSize);
                
                text((size(tmp,1)-3)/2, 0.000001, sprintf('n=%i', sum(any(~isnan(tmp)))))
                
            end
            
            subplot(2,size(DATA2PLOT,2)/2,1)
            t=ylabel(sprintf('Left ROIs\n\nWeight'));
            set(t,'fontsize',FontSize);
            
            subplot(2,size(DATA2PLOT,2)/2,1+size(DATA2PLOT,2)/2)
            t=ylabel(sprintf('Right ROIs\n\nWeight'));
            set(t,'fontsize',FontSize);
            
            print(gcf, fullfile(StartDirectory,  'ROI_Analysis', 'Figures', ...
                [FigName '.tif']), '-dtiff')
            
            
            
            if iRestrict==1
                FigName = ['SVM_Weight_No_Subject_' num2str(size(DATA2PLOT,1)-1) '_Layers_' Analysis(iSVM).name];
            else
                FigName = ['SVM_Weight_Restrict_No_Subject_' num2str(size(DATA2PLOT,1)-1) '_Layers_' Analysis(iSVM).name];
            end
            
            figure('name', FigName,...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible)
            
            
            for iSubROI = 1:size(DATA2PLOT,2)
                
                tmp = squeeze(DATA2PLOT(:,iSubROI,:));
                
                tmp(tmp==0)=NaN;
                
                tmp=fliplr(tmp);
                
                %%
                subplot(2,size(DATA2PLOT,2)/2,iSubROI)
                
                hold on
                grid on
                
                errorbar(0.9, nanmean(tmp(:)), nanstd(tmp(:)), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 3)
                
                shadedErrorBar(2:size(tmp,1), nanmean(tmp(1:end-1,:),2),nansem(tmp(1:end-1,:),2), ...
                    {'LineWidth', 3, 'Color', 'b'}, 1)
                
                tmp2 = nanmean(tmp(1:end-1,:),2)+nansem(tmp(1:end-1,:),2);
                
                axis([0.5 size(tmp,1) 0 max(tmp2(:))])
                
                clear tmp2
                
                set(gca,'tickdir', 'out', ...
                    'xtick', [2 3 size(tmp,1)], ...
                    'xticklabel', ['O';'W';'P'],'fontsize',FontSize)
                
                t=xlabel(SubROIs{iSubROI}(1:end-1));
                set(t,'fontsize',FontSize);
                
                text((size(tmp,1)-3)/2, 0.000001, sprintf('n=%i', sum(any(~isnan(tmp)))))
                
            end
            
            subplot(2,size(DATA2PLOT,2)/2,1)
            t=ylabel(sprintf('Left ROIs\n\nWeight'));
            set(t,'fontsize',FontSize);
            
            subplot(2,size(DATA2PLOT,2)/2,1+size(DATA2PLOT,2)/2)
            t=ylabel(sprintf('Right ROIs\n\nWeight'));
            set(t,'fontsize',FontSize);
            
            print(gcf, fullfile(StartDirectory,  'ROI_Analysis', 'Figures', ...
                [FigName '.tif']), '-dtiff')
            
            
            
            
        end
        
        close all
        
    end
    
end
