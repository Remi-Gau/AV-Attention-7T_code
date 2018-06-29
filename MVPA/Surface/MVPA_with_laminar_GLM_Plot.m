clear
clc

StartDirectory = fullfile(pwd, '..', '..', '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))
cd(StartDirectory)
Get_dependencies('D:\Dropbox')

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

ROIs = {...
    'A1'; ...
    'PT'; ...
    'V1'; ...
    'V2-3';...
    };

%     Analysis(1) = struct('name', 'A Stim VS V Stim', 'class', [1 2]);
%     Analysis(end+1) = struct('name', 'A Stim VS AV Stim', 'class', [1 3]);
%     Analysis(end+1) = struct('name', 'V Stim VS AV Stim', 'class', [2 3]);

Analysis(1) = struct('name', 'A Att VS V Att');
%     Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)', 'class', [6 7]);
%     Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)', 'class', [8 9]);
%     Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)', 'class', [10 11]);

Analysis(end+1) = struct('name', '[Att_A-Att_V]_A VS [Att_A-Att_V]_V');
Analysis(end+1) = struct('name', '[AV-A]_Att_A VS [AV-A]_Att_V');
Analysis(end+1) = struct('name', '[AV-V]_Att_A VS [AV-V]_Att_V');

Analysis(end+1) = struct('name', '[A Att_A VS A Att_V] gen to [V Att_A VS V Att_V]');
Analysis(end+1) = struct('name', '[AV VS A]_Att_A gen to [AV VS A]_Att_V');
Analysis(end+1) = struct('name', '[AV VS V]_Att_A gen to [AV VS V]_Att_V');


NbLayers = 6;

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

SaveSufix = '_results_surf';
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
SaveSufix = [SaveSufix '_FWHM_0.mat'];


%%
for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n\nAnalysing subject %s\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory,...
        'Subjects_Data', ['Subject_' SubjID]);
    
    SaveDir = fullfile(SubjectFolder,'Results', 'SVC');
    
    for iSVM=1:numel(Analysis)
        
        fprintf(' SVM: %s.\n', Analysis(iSVM).name)
        
        for iROI=1:size(ROIs,1)
            
            fprintf('  ROI: %s\n', ROIs{iROI,1})
            
            clear CV Results temp
            
            load(fullfile(SaveDir, ...
                ['SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} SaveSufix]), 'Results');
            
            CV = Results.session(end).rand.perm.CV;
            
            Class_Acc(SubjInd,iSVM,iROI) = nanmean([CV(:).acc]);
            Class_Cst_Acc(SubjInd,iSVM,iROI) = nanmean([CV(:).acc_cst]);
            Class_Cst_Lin(SubjInd,iSVM,iROI) = nanmean([CV(:).acc_lin]);
            
            for iLayer= 1:NbLayers
                Label = [];
                Pred = [];
                for iCV=1:numel(CV)
                    Label(:,end+1) = CV(iCV).layers.results{1}{iLayer}.label;
                    Pred(:,end+1) = CV(iCV).layers.results{1}{iLayer}.pred(:,iLayer);
                end
                MVPA_SubjectsData(SubjInd,iSVM,iROI,iLayer) = mean(mean(Pred==Label,2));
                Acc(iLayer,:) = mean(Pred==Label,2)';
            end
            
            X=repmat(DesMat,size(Acc,2),1);
            Acc=flipud(Acc(:)-.5);
            [B,~,~] = glmfit(X, Acc, 'normal', 'constant', 'off');
            SubjectsBetas(SubjInd,iSVM,iROI,1:size(X,2)) = B;
            
            clear Acc Pred Label B X
            
            clear iLayers
            
        end
        clear iROI
        
    end
    clear iSVM
    
end


% Compile p values and betas
for iSVM = 1:numel(Analysis)
    for iROI=1:numel(ROIs)
        tmp = squeeze(SubjectsBetas(:,iSVM,iROI,:));
        [~,P] = ttest(tmp);
        All_P(iROI,:,iSVM) = P;
        
        MVPA_SubjectsBetas(:,iROI,1:size(tmp,2),iSVM) = tmp;
    end
end
clear P SubjectsBetas




%% Plot
close all
clear DATA

Legends1 = {'', 'Constant', '','','', '', '', '', '', 'Linear'};
Legends2 = {...
    'ROI', ...
    'mean', '(','STD',')', 't value','p value', 'effect size', '', ...
    'mean', '(','STD',')', 't value','p value', 'effect size'};

SavedTxt = fullfile(StartDirectory,'Figures','Generalize_Attention.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, 'MVPA profile\n');
for i=1:length(Legends1)
    fprintf (fid, '%s,', Legends1{i});
end
fprintf (fid, '\n');
for i=1:length(Legends2)
    fprintf (fid, '%s,', Legends2{i});
end
fprintf (fid, '\n');

FigureFolder = fullfile(StartDirectory, 'Figures', strcat(num2str(NbLayers+2), '_layers'));

DATA.WithSubj = 0;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.MVPA = 1;
DATA.Color =  [0 0 0];

Fontsize=12;
Switch = 1;
FigDim = [50, 50, 1300, 700];

for iSVM = 1:numel(Analysis)
    
    figure('position', FigDim, 'name', Analysis(iSVM).name, 'Color', [1 1 1], 'visible', 'on')
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    fprintf (fid, '%s\n', Analysis(iSVM).name);

    for iROI = 1 :numel(ROIs)

        DATA.MIN = 0.35;
        DATA.MAX = 1;   
        subplot(3,4,iROI)
        PlotRectangle(NbLayers,Fontsize,Switch)
        subplot(3,4,iROI)
        DATA.Name = ROIs{iROI};
        DATA.Data = squeeze(MVPA_SubjectsData(:,iSVM,iROI,:));
        DATA.Betas = squeeze(MVPA_SubjectsBetas(:,iROI,1:2,iSVM));
        DATA.Betas(:,2)=DATA.Betas(:,2)*-1;
        DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
        DATA.YLabel = 'Decoding accuracy';
        PlotProfileAndBetas(DATA)
        title(ROIs{iROI})
        
        DATA.InsetLim = [0.5 0.12; -.5 -.12];
        ax = subplot(3,4,iROI+4);
        axis('off')
        DATA.ax = ax.Position;
        DATA.YLabel = 'S Param. est. [a u]';
        DATA.ToPermute = ToPermute; PlotInsetFinal(DATA)
        
        Print2Table(fid, ROIs, iROI, DATA)
        
        DATA.InsetLim = [1 1; 0.1 0.1];
        ax = subplot(3,4,iROI+8);
        axis('off')
        DATA.ax = ax.Position;
        DATA.Betas = [Class_Cst_Acc(:,iSVM,iROI) Class_Cst_Lin(:,iSVM,iROI)];
        DATA.YLabel = 'Decoding accuracy on S Param. Est.';
        PlotInsetFinalAccuracy(DATA)
%         plot([0 1.55], [0.5 0.5], ':k', 'LineWidth', .5)

    end
    
    mtit(strrep(Analysis(iSVM).name,'_Att','_{Att}'), 'xoff', 0, 'yoff', +0.03, 'fontsize', 16)
    
    print(gcf, fullfile(FigureFolder, strcat(Analysis(iSVM).name, '_', num2str(NbLayers), 'Layers.tif')), '-dtiff')
end

fclose (fid);