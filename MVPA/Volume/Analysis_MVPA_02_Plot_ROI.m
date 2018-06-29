clc; clear; close all;

StartDirectory = fullfile(pwd, '..', '..');

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
Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
% Analysis(end+1) = struct('name', 'AV VS A+V');
Analysis(end+1) = struct('name', 'A Att VS V Att');
% Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
% Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');

ROIs= {...
    'V1'
    'V2-3'

    %     'V4'
    %     'V5'
    %     'V1-2-3'
    
    'TE'
    'STGpost'
    };

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

FigDim = [100 100 1500 1000];

Visible = 'on';

n=2;
m=4;

Transparent = 1;

Fontsize = 12;

Scatter = linspace(0,.4,11);

for NbLayers=[6]
    
    for ImgNorm= [0]
        
        for WithQuad= [0]
            
            DesMat = (1:NbLayers)-mean(1:NbLayers);
            if WithQuad
                DesMat = [DesMat' (DesMat.^2)' ones(NbLayers,1)];
            else
                DesMat = [DesMat' ones(NbLayers,1)];
            end
            
            
            
            for WithPerm = [0]
                
                opt.scaling.img.eucledian = 0;
                opt.scaling.img.zscore = ImgNorm;
                opt.scaling.feat.mean = 1;
                opt.scaling.feat.range = 0;
                opt.scaling.feat.sessmean = 0;
                opt.scaling.idpdt = 1;
                
                
                FigureFolder = fullfile(StartDirectory, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));
                cd(FigureFolder)
%                 load(strcat('Data_Block_', num2str(NbLayers), '_Layers', '.mat'), 'AllSubjects_Data')
%                 Include =[];
%                 
%                 AllSubjects_Data(2).name='STGpost';
%                 
%                 for iROI=1:numel(ROIs)
%                     if any(strcmp(ROIs{iROI},{AllSubjects_Data.name}'))
%                         temp = find(strcmp(ROIs{iROI},{AllSubjects_Data.name}'));
%                         Include(:,end+1) = AllSubjects_Data(temp).Include; %#ok<*SAGROW>
%                     end
%                 end
%                 
%                 Include([4 11],:) = [];
                
                Include = true(11,4);
                
                
                for SubjInd = 1:size(SubjectList,1)
                    
                    SubjID = SubjectList(SubjInd,:);
                    
                    fprintf('\n\nAnalysing subject %s\n', SubjID)
                    
                    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
                    
                    SaveDir = fullfile(SubjectFolder, 'Transfer', 'SVM');
                    
                    
                    %%
                    for iSVM=1:numel(Analysis)
                        
                        fprintf('\n SVM: %s.\n', Analysis(iSVM).name)
                        
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
                            
                            SaveSufix = [SaveSufix '_FWHM_6_Layers_' num2str(NbLayers) '.mat'];
                            
                            Save_vol = [Save_vol SaveSufix];
                            
                            load(fullfile(SaveDir, Save_vol), 'Class_Acc', 'Results');
                            AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1)) = squeeze(Class_Acc.TotAcc(end,:,:,:));
                            
                            for iLayer = 2:NbLayers+1
                                Label = cat(2,Results.session(end).rand.perm.CV(:,iLayer).label);
                                Pred = cat(2,Results.session(end).rand.perm.CV(:,iLayer).pred);
                                Acc(iLayer-1,:) = mean(Pred==Label,2)';
                            end
                            
                            X=repmat(DesMat,size(Acc,2),1);
                            Acc=flipud(Acc(:)-.5);
                            [B,~,~] = glmfit(X, Acc, 'normal', 'constant', 'off');
                            
                            AllSubjectsBetas(SubjInd,iSVM,iROI,1:numel(B)) = B;
                            
                            clear Acc Pred Label B X
                        end
                        
                    end
                    
                end
                
                
                %%
                for iSVM = 1:numel(Analysis)
                    for iROI=1:numel(ROIs)
                        tmp = squeeze(AllSubjectsBetas(logical(Include(:,iROI)),iSVM,iROI,:));
                        [~,P] = ttest(tmp);
                        All_P(iROI,:,iSVM) = P;
                        
                        All_Betas(:,iROI,1:size(tmp,2),iSVM) = tmp;
                    end
                end
                clear P
                
                
                %%
                HolmsThresholds = .05*ones(1,numel(ROIs))./(numel(ROIs)+1-(1:numel(ROIs)));
                for iSVM = 1:numel(Analysis)
                    for iP=1:size(All_P,2)
                        tmp = All_P(:,iP,iSVM);
                        [~,I] = ismember(tmp,sort(tmp,'ascend'));
                        Thresholds(:,iP,iSVM) = HolmsThresholds(I);
                    end
                end
                clear All_P
                
                
                %% Permutations
                % Get all possible permutations
                for iSubj=1:size(All_Betas,1)
                    sets{iSubj} = [0 1];
                end
                [a,b,c,d,e,f,g,h,i,j,k] = ndgrid(sets{:});
                cartProd = logical([a(:),b(:),c(:),d(:),e(:),f(:),g(:),h(:),i(:),j(:),k(:)]);
                
                % Do the permutations
                NbPerm=size(cartProd,1);
                NullDist = nan(NbPerm,size(All_Betas,3));
                
                for iSVM = 1:numel(Analysis)
                    for iPerm=1:NbPerm
                        tmp = All_Betas(:,:,:,iSVM);
                        tmp(cartProd(iPerm,:),:,:) = tmp(cartProd(iPerm,:),:,:)*-1;
                        NullDist(iPerm,:) = squeeze(max(abs(mean(tmp)),[],2));
                    end
                    AllNullDist(:,:,iSVM) = sort(NullDist);
                end
                
                
                %% AllSubjects
                
                SaveSufix = '';
                
                if opt.scaling.idpdt
                    SaveSufix = [SaveSufix '-Idpdt']; %#ok<*AGROW>
                end
                if opt.scaling.img.zscore
                    SaveSufix = [SaveSufix '-IMG:ZScore'];
                end
                if opt.scaling.img.eucledian
                    SaveSufix = [SaveSufix '-IMG:Eucl'];
                end
                if opt.scaling.img.zscore==0 && opt.scaling.img.eucledian==0
                    SaveSufix = [SaveSufix '-IMG:None'];
                end
                if opt.scaling.feat.mean
                    SaveSufix = [SaveSufix '-FEAT:MeanCent'];
                end
                if opt.scaling.feat.range
                    SaveSufix = [SaveSufix '-FEAT:Range'];
                end
                if opt.scaling.feat.sessmean
                    SaveSufix = [SaveSufix '-FEAT:SessMeanCent'];
                end
                
                if WithQuad
                    SaveSufix = [SaveSufix '-GLM:Quad'];
                else
                    SaveSufix = [SaveSufix '-GLM:NoQuad'];
                end
                
                if WithPerm
                    SaveSufix = [SaveSufix '-SIG:Perm'];
                else
                    SaveSufix = [SaveSufix '-SIG:NoPerm'];
                end
                
                
                SaveSufix = [SaveSufix '-FWHM6'];
                
                for iROI=1%:numel(ROIs)
                    
                    figure('position', FigDim, 'name', Analysis(iSVM).name, 'Color', [1 1 1], 'visible', Visible)
                    
                    set(gca,'units','centimeters')
                    pos = get(gca,'Position');
                    ti = get(gca,'TightInset');
                    
                    set(gcf, 'PaperUnits','centimeters');
                    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                    set(gcf, 'PaperPositionMode', 'manual');
                    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                    
                    iSubPlot=5;
                    
                    for iSVM = 1:numel(Analysis)
                        
                        subplot(n,m,iSubPlot)
                        
                        DATA.Name = Analysis(iSVM).name;
                        
                        DATA.Data = fliplr(squeeze(AllSubjects(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1))));
                        DATA.Betas = squeeze(AllSubjectsBetas(:,iSVM,iROI,:));
                        
                        
                        DATA.WithPerm = WithPerm;
                        if WithPerm
                            DATA.AllNullDist = AllNullDist(:,:,iSVM);
                        else
                            DATA.Thresholds = squeeze(Thresholds(iROI,:,iSVM));
                        end
                        
                        
                        if n==2 && m==3
                            switch iSubPlot
                                case 1
                                    XY = [0.165 0.605 .04 .07];
                                case 2
                                    XY = [0.445 0.605 .04 .07];
                                case 3
                                    XY = [0.73 0.605 .04 .07];
                                case 4
                                    XY = [0.165 0.135 .04 .07];
                                case 5
                                    XY = [0.445 0.135 .04 .07];
                            end
                            
                        elseif n==2 && m==2
                            switch iSubPlot
                                case 1
                                    XY = [0.2 0.61 .04 .07];
                                case 2
                                    XY = [0.65 0.61 .04 .07];
                                case 3
                                    XY = [0.2 0.14 .04 .07];
                                case 4
                                    XY = [0.65 0.14 .04 .07];
                            end
                            
                        elseif n==2 && m==4
                            switch iSubPlot
                                case 1
                                    XY = [0.2 0.61 .04 .07];
                                case 2
                                    XY = [0.4 0.61 .04 .07];
                                case 3
                                    XY = [0.6 0.61 .04 .07];
                                case 4
                                    XY = [0.8 0.61 .04 .07];                                    
                                case 5
                                    XY = [0.175 0.15 .04 .07];
                                case 6
                                    XY = [0.385 0.15 .04 .07];
                                case 7
                                    XY = [0.59 0.15 .04 .07];
                                case 8
                                    XY = [0.8 0.15 .04 .07];                                    
                            end
                        end
                        
                        DATA.XY= XY;
                        DATA.Scatter = Scatter;
                        
                        DATA.MVPA = 1;
                        
                        PlotProfileAndBetas(DATA)
                        
                        
                        iSubPlot = iSubPlot + 1;
                        
                    end
                    
                    SaveSufix = strrep(SaveSufix, '_', '-');
                    
                    t = mtit([ROIs{iROI} SaveSufix], 'xoff', 0, 'yoff', +0.05, ...
                        'fontsize',Fontsize+4);
                    
                    Name = [ROIs{iROI} SaveSufix];
                    
                    print(gcf, fullfile(StartDirectory, 'Figures', 'MVPA', 'vol', ...
                        [Name '_' num2str(NbLayers) 'Layers.pdf']), '-dpdf')
                    
                    %     print(gcf, fullfile(StartDirectory, 'Figures', 'MVPA', 'vol', ...
                    %         [Name '_' num2str(NbLayers) 'Layers.tif']), '-dtiff')
                    
%                     close all
                    
                end
                
               % clear AllNullDist AllSubjectsBetas
            end
        end
    end
end