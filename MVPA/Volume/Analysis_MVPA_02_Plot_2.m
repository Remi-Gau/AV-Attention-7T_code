clc; clear; close all;

StartDirectory = pwd;

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


FigDim = [100 100 1500 1000];

Visible = 'off';

Analysis(1) = struct('name', 'A Stim VS V Stim');
Analysis(end+1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');

ROIs= {...
    'V1'
    'V2'
    'V3'
    'V4'
    'V5'
    'V1-2-3'
    
    'TE1.0'
    'TE1.1'
    'TE1.2'
    'TE'
    'STGpost'
    };

Sel = {'0','3', '6'};

NbLayers = 6;

FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(6), '_layers'));
cd(FigureFolder)
load(strcat('Data_Block_', num2str(6), '_Layers', '.mat'), 'AllSubjects_Data')
Include =[];

AllSubjects_Data(10).name='STGpost';

for iROI=1:numel(ROIs)
    if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
        temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
        Include(:,end+1) = AllSubjects_Data(temp).Include;
    end
end

Include([4 11],:) = [];



for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n\nAnalysing subject %s\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    SaveDir = fullfile(SubjectFolder, 'Transfer', 'SVM');
    
    
    %%
    
    for iSVM=1:numel(Analysis)
        
        fprintf('\n SVM: %s.\n', Analysis(iSVM).name)
        
        for iROI=1:numel(ROIs)
            
            Norm = 6;
            
            switch Norm
                case 1
                    opt.scaling.img.eucledian = 1;
                    opt.scaling.img.zscore = 0;
                    opt.scaling.feat.mean = 0;
                    opt.scaling.feat.range = 0;
                    opt.scaling.feat.sessmean = 1;
                case 2
                    opt.scaling.img.eucledian = 1;
                    opt.scaling.img.zscore = 0;
                    opt.scaling.feat.mean = 0;
                    opt.scaling.feat.range = 1;
                    opt.scaling.feat.sessmean = 0;
                case 3
                    opt.scaling.img.eucledian = 1;
                    opt.scaling.img.zscore = 0;
                    opt.scaling.feat.mean = 1;
                    opt.scaling.feat.range = 0;
                    opt.scaling.feat.sessmean = 0;
                case 4
                    opt.scaling.img.eucledian = 0;
                    opt.scaling.img.zscore = 1;
                    opt.scaling.feat.mean = 0;
                    opt.scaling.feat.range = 0;
                    opt.scaling.feat.sessmean = 1;
                case 5
                    opt.scaling.img.eucledian = 0;
                    opt.scaling.img.zscore = 1;
                    opt.scaling.feat.mean = 0;
                    opt.scaling.feat.range = 1;
                    opt.scaling.feat.sessmean = 0;
                case 6
                    opt.scaling.img.eucledian = 0;
                    opt.scaling.img.zscore = 1;
                    opt.scaling.feat.mean = 1;
                    opt.scaling.feat.range = 0;
                    opt.scaling.feat.sessmean = 0;
            end
            
            Idpdt = 1;
            opt.scaling.idpdt = Idpdt;
            
            for iSel = 2:length(Sel)
                
                Save_vol = [...
                    'SVM_' Analysis(iSVM).name...
                    '_ROI_' ROIs{iROI}];
                
                SaveSufix = '_results_vol';
                SaveSufix2 = '_results_vol_FixedC';
                
                if opt.scaling.idpdt
                    SaveSufix = [SaveSufix '_Idpdt'];
                    SaveSufix2 = [SaveSufix2 '_Idpdt'];
                end
                if opt.scaling.img.zscore
                    SaveSufix = [SaveSufix '_ZScore'];
                    SaveSufix2 = [SaveSufix2 '_ZScore'];
                end
                if opt.scaling.img.eucledian
                    SaveSufix = [SaveSufix '_Eucl'];
                    SaveSufix2 = [SaveSufix2 '_Eucl'];
                end
                if opt.scaling.feat.mean
                    SaveSufix = [SaveSufix '_MeanCent'];
                    SaveSufix2 = [SaveSufix2 '_MeanCent'];
                end
                if opt.scaling.feat.range
                    SaveSufix = [SaveSufix '_Range'];
                    SaveSufix2 = [SaveSufix2 '_Range'];
                end
                if opt.scaling.feat.sessmean
                    SaveSufix = [SaveSufix '_SessMeanCent'];
                    SaveSufix2 = [SaveSufix2 '_SessMeanCent'];
                end
                
                SaveSufix = [SaveSufix '_FWHM_' Sel{iSel} '_Layers_' NbLayers '.mat'];
                SaveSufix2 = [SaveSufix2 '_FWHM_' Sel{iSel} '_Layers_' num2str(NbLayers) '.mat'];
                
                Save_vol2 = [Save_vol SaveSufix2];
                Save_vol = [Save_vol SaveSufix];
                
                
%                 load(fullfile(SaveDir, Save_vol), 'Class_Acc');
%                 AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1),iSel,1) = squeeze(Class_Acc.TotAcc(end,:,:,:)); %#ok<SAGROW>
                
                load(fullfile(SaveDir, Save_vol2), 'Class_Acc');
                AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1),iSel,1) = squeeze(Class_Acc.TotAcc(end,:,:,:)); %#ok<SAGROW>
                
            end
            
        end
        
    end
    
end



%% AllSubjects

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

mn = length(ROIs);
n  = round(mn^0.4);
m  = ceil(mn/n);

Transparent = 1;

Norm =  6;

switch Norm
    case 1
        opt.scaling.img.eucledian = 1;
        opt.scaling.img.zscore = 0;
        opt.scaling.feat.mean = 0;
        opt.scaling.feat.range = 0;
        opt.scaling.feat.sessmean = 1;
    case 2
        opt.scaling.img.eucledian = 1;
        opt.scaling.img.zscore = 0;
        opt.scaling.feat.mean = 0;
        opt.scaling.feat.range = 1;
        opt.scaling.feat.sessmean = 0;
    case 3
        opt.scaling.img.eucledian = 1;
        opt.scaling.img.zscore = 0;
        opt.scaling.feat.mean = 1;
        opt.scaling.feat.range = 0;
        opt.scaling.feat.sessmean = 0;
    case 4
        opt.scaling.img.eucledian = 0;
        opt.scaling.img.zscore = 1;
        opt.scaling.feat.mean = 0;
        opt.scaling.feat.range = 0;
        opt.scaling.feat.sessmean = 1;
    case 5
        opt.scaling.img.eucledian = 0;
        opt.scaling.img.zscore = 1;
        opt.scaling.feat.mean = 0;
        opt.scaling.feat.range = 1;
        opt.scaling.feat.sessmean = 0;
    case 6
        opt.scaling.img.eucledian = 0;
        opt.scaling.img.zscore = 1;
        opt.scaling.feat.mean = 1;
        opt.scaling.feat.range = 0;
        opt.scaling.feat.sessmean = 0;
end


Idpdt = 1;

opt.scaling.idpdt=Idpdt;

for iSel = 2:length(Sel)
    
    SaveSufix = '_FixedC';
    
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
    
    for iSVM = 1:numel(Analysis)
        
        
        figure('position', FigDim, 'name', Analysis(iSVM).name, 'visible', Visible)
        
        set(gca,'units','centimeters')
        pos = get(gca,'Position');
        ti = get(gca,'TightInset');
        
        set(gcf, 'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        
        iSubPlot=1;
        
        for iROI=1:numel(ROIs)
            
            subplot(n,m,iSubPlot)
            
            hold on
            grid on

            Data = squeeze(AllSubjects(logical(Include(:,iROI)),iSVM,iROI,1:(NbLayers+1),iSel,1));
            
            for SubjInd = 1:size(Data,1)
                plot(1:NbLayers, fliplr(Data(SubjInd,2:(NbLayers+1))), '--', ...
                    'LineWidth', .5, 'Color', COLOR_Subject(SubjInd,:), ...
                    'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:));
            end
            clear SubjInd
            
            shadedErrorBar(1:NbLayers, nanmean(fliplr(Data(:,2:(NbLayers+1)))),nansem(fliplr(Data(:,2:(NbLayers+1)))), ...
                {'LineWidth', 3, 'Color', [0 0 1]}, Transparent)
            
            plot([1 numel(ROIs)], [0.5 0.5], '--k', 'linewidth', 2)
            
            axis([0.5 NbLayers+.5 0 1])
            
            set(gca,'tickdir', 'out', ...
                'xtick', 1:NbLayers, ...
                'xticklabel', 1:NbLayers, ...
                'ytick', 0:.1:1, ...
                'yticklabel', 0:.1:1,...
                'fontsize', 8)
            
            title(ROIs{iROI},'fontsize', 10)
            
            %     legend(char({...
            %         'Res: 2.5 mm - 3 mm FWHM - SessMeanCt';...
            %         'Res: 2.5 mm - 3 mm FWHM - Range';...
            %         'Res: 2.5 mm - 3 mm FWHM - MeanCt';...
            %         }), 'Location', 'South')
            
            iSubPlot = iSubPlot + 1;
            
        end
        
        SaveSufix = strrep(SaveSufix, '_', '-');
        
        mtit([Analysis(iSVM).name '-FWHM ' Sel{iSel} SaveSufix], 'xoff', 0, 'yoff', +0.02)
        
        Name = [strrep(Analysis(iSVM).name, ' ', '_') '_FWHM_' Sel{iSel} SaveSufix '.pdf'];
        Name = strrep(Name, '(', '-');
        Name = strrep(Name, ')', '');
        print(gcf, fullfile(StartDirectory, 'Figures', 'MVPA', 'vol', ...
            Name), '-dpdf')
        
        close all
        
    end
    
end


%% 

cd(fullfile(StartDirectory, 'Figures', 'MVPA', 'vol'))

close all

Command = [];


for iSVM = 1:numel(Analysis)
    
    Norm =  6;
    
    switch Norm
        case 1
            opt.scaling.img.eucledian = 1;
            opt.scaling.img.zscore = 0;
            opt.scaling.feat.mean = 0;
            opt.scaling.feat.range = 0;
            opt.scaling.feat.sessmean = 1;
        case 2
            opt.scaling.img.eucledian = 1;
            opt.scaling.img.zscore = 0;
            opt.scaling.feat.mean = 0;
            opt.scaling.feat.range = 1;
            opt.scaling.feat.sessmean = 0;
        case 3
            opt.scaling.img.eucledian = 1;
            opt.scaling.img.zscore = 0;
            opt.scaling.feat.mean = 1;
            opt.scaling.feat.range = 0;
            opt.scaling.feat.sessmean = 0;
        case 4
            opt.scaling.img.eucledian = 0;
            opt.scaling.img.zscore = 1;
            opt.scaling.feat.mean = 0;
            opt.scaling.feat.range = 0;
            opt.scaling.feat.sessmean = 1;
        case 5
            opt.scaling.img.eucledian = 0;
            opt.scaling.img.zscore = 1;
            opt.scaling.feat.mean = 0;
            opt.scaling.feat.range = 1;
            opt.scaling.feat.sessmean = 0;
        case 6
            opt.scaling.img.eucledian = 0;
            opt.scaling.img.zscore = 1;
            opt.scaling.feat.mean = 1;
            opt.scaling.feat.range = 0;
            opt.scaling.feat.sessmean = 0;
    end
    
    
    Idpdt = 1;
    
    opt.scaling.idpdt=Idpdt;
    
    for iSel = 2:length(Sel)
        
        for i=1
            
            if i==1
                SaveSufix = '_FixedC';
            else
                SaveSufix = [];
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
            
            SaveSufix = strrep(SaveSufix, '_', '-');
            
            Name = [strrep(Analysis(iSVM).name, ' ', '_') '_FWHM_' Sel{iSel} SaveSufix '.pdf'];
            Name = strrep(Name, '(', '-');
            Name = strrep(Name, ')', '');
            
            Command = [Command ' ' Name];
        end
    end
end


system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, ['MVPA_vol_ALL.pdf']) ' ' Command])



