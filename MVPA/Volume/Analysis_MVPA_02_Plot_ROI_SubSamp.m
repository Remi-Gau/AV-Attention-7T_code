clc; clear; close all;

StartDirectory = fullfile(pwd, '..', '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))

FigDir = '/home/rxg243/Dropbox/PhD/Experiments/AV_Integration_7T/Figures/MVPA/vol/subsamp';
[~,~,~] = mkdir(FigDir);

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

Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');

Analysis(end+1) = struct('name', 'A Att VS V Att');


ROIs= {...
    'V1'
    'V2-3'
    
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
m=3;

Transparent = 1;

Fontsize = 12;

Scatter = linspace(0,.4,11);

for NbLayers=[6]
    
    for ImgNorm= [0 1]
        
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
                            
                            SaveSufix = '_results_vol_FixedC_SubSamp';
                            
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
                            clear tmp
                            for i=1:size(Class_Acc,1)
                                for j=1:size(Class_Acc,2)
                                    tmp(:,i,j) = squeeze(Class_Acc(i,j).TotAcc(end,:,:,:));
                                end
                            end
                            
                            AllSubjects(SubjInd,iSVM,iROI,1:(NbLayers+1),1:size(Class_Acc,2))= mean(tmp,2); clear tmp
                            
                            for i=1:size(Class_Acc,1)
                                for j=1:size(Class_Acc,2)
                                    for iLayer = 2:NbLayers+1
                                        Label = cat(2,Results.session(end).rand.perm.SubSamp{i,j}.CV(:,iLayer).label);
                                        Pred = cat(2,Results.session(end).rand.perm.SubSamp{i,j}.CV(:,iLayer).pred);
                                        Acc(iLayer-1,:,j,i) = mean(Pred==Label,2)';
                                    end
                                end
                            end
                            
                            Acc = mean(Acc,4);
                            
                            for j=1:size(Class_Acc,2)
                                X=repmat(DesMat,size(Acc,2),1);
                                tmp = Acc(:,:,j);
                                tmp=flipud(tmp(:)-.5);
                                [B,~,~] = glmfit(X, tmp, 'normal', 'constant', 'off');
                                
                                AllSubjectsBetas(SubjInd,iSVM,iROI,1:numel(B),j) = B;
                                
                                clear tmp Pred Label B X
                            end
                            clear Acc
                            
                            
                            
                        end
                        
                    end
                    
                end
                
                
                %%
                for iSVM = 1:numel(Analysis)
                    for iROI=1:numel(ROIs)
                        for i=1:size(AllSubjectsBetas,5)
                            
                            tmp = squeeze(AllSubjectsBetas(logical(Include(:,iROI)),iSVM,iROI,:,i));
                            [~,P] = ttest(tmp);
                            All_P(iROI,:,iSVM,i) = P; %#ok<*SAGROW>
                            
                            All_Betas(:,iROI,1:size(tmp,2),iSVM,i) = tmp;
                            
                        end
                    end
                end
                clear P
                
                
                
                %%
                HolmsThresholds = .05*ones(1,numel(ROIs))./(numel(ROIs)+1-(1:numel(ROIs)));
                for iSVM = 1:numel(Analysis)
                    for iP=1:size(All_P,2)
                        for i=1:size(All_P,4)
                            tmp = All_P(:,iP,iSVM,i);
                            [~,I] = ismember(tmp,sort(tmp,'ascend'));
                            Thresholds(:,iP,iSVM,i) = HolmsThresholds(I);
                        end
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
                    for i=1:size(All_Betas,5)
                        for iPerm=1:NbPerm
                            tmp = All_Betas(:,:,:,iSVM,i);
                            tmp(cartProd(iPerm,:),:,:) = tmp(cartProd(iPerm,:),:,:)*-1;
                            NullDist(iPerm,:) = squeeze(max(abs(mean(tmp)),[],2));
                        end
                        AllNullDist(:,:,iSVM,i) = sort(NullDist);
                    end
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
                
                
                SaveSufix = [SaveSufix '-FWHM6_SubSamp'];
                
                for iROI=1:numel(ROIs)
                    
                    figure('position', FigDim, 'name', Analysis(iSVM).name, 'Color', [1 1 1], 'visible', Visible)
                    
                    set(gca,'units','centimeters')
                    pos = get(gca,'Position');
                    ti = get(gca,'TightInset');
                    
                    set(gcf, 'PaperUnits','centimeters');
                    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                    set(gcf, 'PaperPositionMode', 'manual');
                    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                    
                    SubPlot = [1 4 2 5 3 6];
                    
                    iSubPlot=1;
                    
                    for iSVM = 1:numel(Analysis)
                        
                        for iSubSamp = 1:2
                            
                            subplot(n,m,SubPlot(iSubPlot))
                            
                            DATA.Name = Analysis(iSVM).name;
                            
                            DATA.Data = fliplr(squeeze(AllSubjects(logical(Include(:,iROI)),iSVM,iROI,2:(NbLayers+1),iSubSamp)));
                            DATA.Betas = squeeze(AllSubjectsBetas(:,iSVM,iROI,:,iSubSamp));
                            
                            DATA.WithPerm = WithPerm;
                            if WithPerm
                                DATA.AllNullDist = AllNullDist(:,:,iSVM,iSubSamp);
                            else
                                DATA.Thresholds = squeeze(Thresholds(iROI,:,iSVM,iSubSamp));
                            end
                            
                            switch SubPlot(iSubPlot)
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
                                case 6
                                    XY = [0.73 0.135 .04 .07];
                            end
                            
                            DATA.XY= XY;
                            DATA.Scatter = Scatter;
                            
                            DATA.WithSubj = 1;
                            
                            DATA.MVPA = 1;
                            
                            PlotProfileAndBetas(DATA)
                            
                            iSubPlot = iSubPlot + 1;
                            
                        end
                    end
                    
                    SaveSufix = strrep(SaveSufix, '_', '-');
                    
                    
                    t = mtit([ROIs{iROI} SaveSufix], 'xoff', 0, 'yoff', +0.05, ...
                        'fontsize',Fontsize+4);
                    
                    SaveSufix = strrep(SaveSufix, ':', '-');
                    Name = [ROIs{iROI} SaveSufix];
                    
                    print(gcf, fullfile(FigDir, [Name '_' num2str(NbLayers) 'Layers.pdf']), '-dpdf')
                    
                    %     print(gcf, fullfile(StartDirectory, 'Figures', 'MVPA', 'vol', ...
                    %         [Name '_' num2str(NbLayers) 'Layers.tif']), '-dtiff')
                    
                    %                     close all
                    
                end
                
                % clear AllNullDist AllSubjectsBetas
            end
        end
    end
end