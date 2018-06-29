clc; clear; close all;

StartDirectory = pwd;

addpath(genpath(fullfile(StartDirectory, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    %'06';...
    %'07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    %'14';...
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
    
    'TE1.0'
    'TE1.1'
    'TE1.2'
    'TE'
    };

Sel = {'0','3', '6'};

NbLayers = 6;


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n\nAnalysing subject %s\n', SubjID)
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    SaveDir = fullfile(SubjectFolder, 'Transfer', 'SVM');
    
    
    %%
    
    for iSVM=1:numel(Analysis)
        
        fprintf('\n SVM: %s.\n', Analysis(iSVM).name)
        
        for iROI=1:numel(ROIs)
            
            for Norm =  [1:6]
                
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
                
                for FeatureRestriction = 1
                    
                    switch FeatureRestriction
                        case 1
                            opt.fs.do = 0;
                            opt.rfe.do = 0;
                        case 2
                            opt.fs.do = 0;
                            opt.rfe.do = 1;
                        case 3
                            opt.fs.do = 1;
                            opt.rfe.do = 0;
                    end
                    
                    for iSel = 1:length(Sel)
                        
                        for Idpdt = 0:1
                            
                            opt.scaling.idpdt = Idpdt;
                            
                            Save_vol = [...
                                'SVM_' Analysis(iSVM).name...
                                '_ROI_' ROIs{iROI}];
                            
                            %                             if Norm==2 && Idpdt==1
                            %                                 SaveSufix = '_results_ds';
                            %                             else
                            SaveSufix = '_results_vol';
                            %                             end
                            
                            if opt.fs.do
                                SaveSufix = [SaveSufix '_FS']; %#ok<*AGROW>
                            end
                            if opt.rfe.do
                                SaveSufix = [SaveSufix '_RFE'];
                            end
                            %                             if opt.permutation.test
                            %                                 SaveSufix = [SaveSufix '_Perm'];
                            %                             end
                            %                             if opt.session.curve
                            %                                 SaveSufix = [SaveSufix '_Lear'];
                            %                             end
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
                            
                            SaveSufix = [SaveSufix '_FWHM_' Sel{iSel} '_Layers_' NbLayers '.mat'];
                            
                            Save_vol = [Save_vol SaveSufix];
                            
                            if exist(fullfile(SaveDir, Save_vol), 'file');
                                load(fullfile(SaveDir, Save_vol), 'Class_Acc');
                                AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1),iSel,Idpdt+1,Norm) = squeeze(Class_Acc.TotAcc(end,:,:,:)); %#ok<SAGROW>
                            else
                                warning('The file %s does not exist.', fullfile(SaveDir, Save_vol))
                                AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1),iSel,Idpdt+1,Norm) = nan(1,NbLayers+1);
                            end
                            
                        end
                        
                    end
                    
                    
                end
                
            end
            
        end
        
    end
    
end

%%
close all

MEAN = squeeze(mean(AllSubjects,1));

Best= [];

for iSVM = 1:size(MEAN,1)
    
    figure('position', FigDim, 'name', Analysis(iSVM).name, 'visible', 'on')
    
    for iROI = 1:size(MEAN,2)

        ToPLot = squeeze(MEAN(iSVM,iROI,2:end,:,:,:));
        ToPLot = reshape(ToPLot, [3,36])';
        [Y,I] = max(mean(ToPLot,2));
        
        Best = [Best I];
        
        subplot(1,size(MEAN,2),iROI)
        hold on
        
        imagesc(ToPLot, [0 1])
        rectangle('Position', [0.5 I-.5 3 1], 'edgecolor', 'r', 'linewidth', 2)
        
        axis([0.5 NbLayers+.5 0.5 36.5])
        
        xlabel('Layers')
        
        if iROI==1
        set(gca,'tickdir', 'out', ...
            'xtick', 1:NbLayers, ...
            'xticklabel', 1:NbLayers, ...
            'ytick', 1:3:36, ...
            'yticklabel', 1:3:36,...
            'fontsize', 8)
        else
            set(gca,'tickdir', 'out', ...
                'xtick', 1:NbLayers, ...
                'xticklabel', 1:NbLayers, ...
                'ytick', [], ...
                'yticklabel', [],...
                'fontsize', 8)
        end
        
        title(ROIs{iROI})
        
    end
    
     mtit(Analysis(iSVM).name, 'xoff', 0, 'yoff', +0.02)
end

tabulate(Best);

Labels = [num2str([1:36]') repmat([' FWHM=0';' FWHM=3';' FWHM=6'],[12,1]) repmat([' Id=0';' Id=1'], [18,1]) ...
    repmat([' IMG=Eucl';' IMG=Eucl';' IMG=Eucl';' IMG=Zsco';' IMG=Zsco';' IMG=Zsco'], [6,1]) ...
    repmat([' FEAT=SeMea';' FEAT=Range';' FEAT=Mean '], [12,1])];

Labels(unique(Best),:);

%%
clc
T = tabulate(Best);
T = [T(:,2) repmat([1:3]', [12,1])];
FWHM = [];
for i=1:size(T,1)
    FWHM = [FWHM ; ones(T(i,1),1)*T(i,2)];
end
tabulate(FWHM)

%%
clc
T = tabulate(Best);
T = [T(:,2) repmat([1:2]', [18,1])];
Idpt = [];
for i=1:size(T,1)
    Idpt = [Idpt ; ones(T(i,1),1)*T(i,2)];
end
tabulate(Idpt)

%%
clc
T = tabulate(Best);
T = [T(:,2) repmat([1;1;1;2;2;2], [6,1])];
IMG = [];
for i=1:size(T,1)
    IMG = [IMG ; ones(T(i,1),1)*T(i,2)];
end
tabulate(IMG)

%%
clc
T = tabulate(Best);
T = [T(:,2) repmat([1;2;3], [12,1])];
FEAT = [];
for i=1:size(T,1)
    FEAT = [FEAT ; ones(T(i,1),1)*T(i,2)];
end
tabulate(FEAT)

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

for Norm =  1:6
    
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
    
    
    for Idpdt = 0:1
        
        opt.scaling.idpdt=Idpdt;
        
        for iSel = 1:length(Sel)
            
            SaveSufix = [];
            
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
                    
                    Data = squeeze(AllSubjects(:,iSVM,iROI,1:(NbLayers+1),iSel,Idpdt+1,Norm));
                    
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
        
    end
end

%% Smoothing

cd(fullfile(StartDirectory, 'Figures', 'MVPA', 'vol'))

close all

Command = [];


for iSVM = 1:numel(Analysis)
    
    for Norm =  1:6
        
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
        
        
        for Idpdt = 0:1
            
            opt.scaling.idpdt=Idpdt;
            
            for iSel = 1:length(Sel)
                
                SaveSufix = [];
                
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
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, ['MVPA_vol_Smoothing.pdf']) ' ' Command])



%% Idependent scaling

cd(fullfile(StartDirectory, 'Figures', 'MVPA', 'vol'))

close all

Command = [];


for iSVM = 1:numel(Analysis)
    
    for Norm =  1:6
        
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
        
        for iSel = 1:length(Sel)
            
            for Idpdt = 0:1
                
                opt.scaling.idpdt=Idpdt;
                
                SaveSufix = [];
                
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
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, ['MVPA_vol_IdpdtScaling.pdf']) ' ' Command])



%% IMG scaling

cd(fullfile(StartDirectory, 'Figures', 'MVPA', 'vol'))

close all

Command = [];


for iSVM = 1:numel(Analysis)
    
    for iSel = 1:length(Sel)
        
        for Idpdt = 0:1
            
            for Norm =  [1 4 2 5 3 6]
                
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
                
                
                
                opt.scaling.idpdt=Idpdt;
                
                SaveSufix = [];
                
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
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, ['MVPA_vol_ImgScaling.pdf']) ' ' Command])



%% IMG scaling

cd(fullfile(StartDirectory, 'Figures', 'MVPA', 'vol'))

close all

Command = [];


for iSVM = 1:numel(Analysis)
    
    for iSel = 1:length(Sel)
        
        for Idpdt = 0:1
            
            for Norm =  1:6
                
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
                
                
                
                opt.scaling.idpdt=Idpdt;
                
                SaveSufix = [];
                
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
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartDirectory, ['MVPA_vol_FeatScaling.pdf']) ' ' Command])

