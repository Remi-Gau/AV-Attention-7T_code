clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces'); %#ok<*NASGU>

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
%     '14';...
    '15';...
    '16'
    ];

% Color map
X = 0:0.001:1;
R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
ColorMap = [R' G' B'];
clear R G B


NbLayers = 7;
NbLayers = NbLayers+1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);


FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

Print = 0;

% ROI_Groups = {1:4;5:8;9:12;13:16;17:20;21:24;25:28};
% ROI_Groups_names = {'A1-PT-V123';'V123-ActDeact';'A1-PT-A-ActDeact';'V123-A-ActDeact';...
%     'A1-PT-V-ActDeact';'V123-V-ActDeact';'A1-PT-ActDeact'};
ROI_Groups = {1:4};
ROI_Groups_names = {'A1-PT-V123'};

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')

%%
for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'FFX_Block');
    
    Data_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'BetaMapping','8Surf');
    
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Surfaces');
    [~,~,~]=mkdir(Results_Folder);
    
    [~,~,~]=mkdir(fullfile(Results_Folder,'Cdtions'));
    [~,~,~]=mkdir(fullfile(Results_Folder,'Att'));
    [~,~,~]=mkdir(fullfile(Results_Folder,'CrossSens'));
    
    %% Creates a cell that lists the names of the beta images as well as their column number in the design matrix
    load(fullfile(GLM_Folder, 'SPM.mat'))
    
    BetaNames = char(SPM.xX.name');
    
    BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
        strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
        2);
    
    BetaOfInterest = find(BetaOfInterest);
    
    BetaNames(:,1:6)=[];
    
    clear SPM
    
    
    %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
    %% Read features
    fprintf(' Reading VTKs\n')
    
    % Format for reading the vertices from the VTK file
    Spec = repmat('%f ', 1, NbLayers);
    
    NbVertices = nan(1,2);
    
    % For the 2 hemispheres
    for hs = 1:2
        
        if hs==1
            fprintf('   Left hemipshere\n')
            HsSufix = 'l';
        else
            fprintf('   Right hemipshere\n')
            HsSufix = 'r';
        end
        
        Betas = dir(fullfile(Data_Folder, ['Beta*' HsSufix 'cr.vtk']));
        
        FeatureSaveFile = fullfile(Data_Folder,[ 'Subj_' SubjID '_features_' HsSufix 'hs_' ...
            num2str(NbLayers) '_surf.mat']);
        
        InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
            ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
        
        [Vertex,~,~] = read_vtk(InfSurfFile, 0, 1);
        NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData; %#ok<*SAGROW>
        else
            
            AllMapping = nan(NbVertices(hs),NbLayers,size(Betas,1));
            
            fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
            
            parfor iBeta = 1:size(Betas,1)
                
                A = fileread(fullfile(Data_Folder,Betas(iBeta).name)); % reads file quickly
                B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
                
                C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
                Mapping = cell2mat(C); %clear C
                
                if size(Mapping,1)~=(NbVertices(hs))  %#ok<*PFBNS>
                    error('A VTK file has wrong number of vertices:\n%s', fullfile(Data_Folder,Betas(iBeta).name))
                end
                
                AllMapping(:,:,iBeta) = Mapping;
                
                fprintf(1,'\b.\n');
                
            end
            fprintf(1,'\b]\n');
            
            A = AllMapping==0;
            A = squeeze(any(A,2));
            A = ~any(A,2);
            VertexWithData = find(A);
            clear A
            
            AllMapping = AllMapping(VertexWithData,:,:);
            
            VertexWithDataHS{hs} = VertexWithData;
            
            save(FeatureSaveFile,'Vertex','Face','AllMapping','VertexWithData', '-v7.3')
        end
        
        clear Betas InfSurfFile Mapping A B C
        
        
        %% Run GLMs for A, V and AV
        for CondInd = 1:length(Conditions_Names)/2 % For each Condition
            
            Beta2Sel = [];
            Beta2Sel2 = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd} ' - Block ' num2str(BlockInd) '*bf(1)']))];
                Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd+3} ' - Block ' num2str(BlockInd) '*bf(1)']))];
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
            
            Features = mean(cat(4,AllMapping(:,2:end-1,Beta2Sel), AllMapping(:,2:end-1,Beta2Sel2)),4); %#ok<*FNDSB>
            FeaturesCdtion{CondInd,hs} = mean(Features,3);
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaCdtion{hs}(:,:,CondInd) = pinv(X)*Y;
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel Beta2Sel2 X Y Mapping
            
        end
        
        
        %% Run GLMs for Attention
        Cond2Contrast = {...
            1, 4;...
            2, 5;...
            3, 6;...
            };
        
        for AttCondInd = 1:size(Cond2Contrast,1) % For each Condition
            
            Beta2Sel = [];
            Beta2Sel2 = [];
            for CondInd=1:numel(Cond2Contrast{AttCondInd,1})
                for BlockInd = 1:3
                    Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                        [Conditions_Names{Cond2Contrast{AttCondInd,1}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
                    Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                        [Conditions_Names{Cond2Contrast{AttCondInd,2}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>.
                end
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
            
            Features = AllMapping(:,2:end-1,Beta2Sel) - ...
                AllMapping(:,2:end-1,Beta2Sel2);
            FeaturesAtt{AttCondInd,hs} = mean(Features,3);
            
            FeaturesAll(:,:,:,AttCondInd) = Features;
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaAtt{hs}(:,:,AttCondInd) = pinv(X)*Y;
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        Features = mean(FeaturesAll,4);
        FeaturesAtt{size(Cond2Contrast,1)+1,hs} = mean(Features,3);
        
        X=repmat(DesMat,size(Features,3),1);
        
        Y = shiftdim(Features,1);
        Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
        
        BetaAtt{hs}(:,:,end+1) = pinv(X)*Y;
        %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
        
        clear FeaturesAll
        
        
        %% Run GLMs for cross-sensory influence
        Cond2Contrast = {...
            [3 6], [1 4];...
            [3 6], [2 5];...
            };
        
        CrossSensCondNames = {...
            'AV-A'; ...
            'AV-V'};
        
        for CrossSensCondInd = 1:size(CrossSensCondNames,1) % For each Condition
            
            for CondInd=1:numel(Cond2Contrast{CrossSensCondInd,1})
                Beta2Sel = [];
                Beta2Sel2 = [];
                for BlockInd = 1:3
                    Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                        [Conditions_Names{Cond2Contrast{CrossSensCondInd,1}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))];
                    Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                        [Conditions_Names{Cond2Contrast{CrossSensCondInd,2}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))];
                end
                Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
                Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
                
                Features(:,:,:,CondInd) = AllMapping(:,2:end-1,Beta2Sel) - ...
                    AllMapping(:,2:end-1,Beta2Sel2);
            end
            
            Features = mean(Features,4);
            
            FeaturesCrossMod{CrossSensCondInd,hs} = mean(Features,3);
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaCrossSens{hs}(:,:,CrossSensCondInd) = pinv(X)*Y;
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        clear CondInd AttCondInd CrossSensCondInd CrossSensCondNames AttCondNames BlockInd CondInd
        
    end
    
    
    if any(NbVertex ~= NbVertices)
        %         NbVertex
        %         NbVertices
        error('The number of vertices does not match.')
    end
    
    close all
    
    
    
    %% Profiles stim = f(Percentile of V)
    close all
    Cdt =[2 1;2 2;2 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
            end
            
            for iCdt = 1:3
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesCdtion{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesCdtion{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp} %1:4 %numel(ROI)
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    All_Profiles_V{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_V{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                        
                        hold on
                        
                        colormap(ColorMap);
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-6 6])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if iROI==ROI_Groups{iRoiGrp}(1)
                            title([Name{Cdt(iCdt,2)} ' stimulus'])
                        end
                        
                        
                        if iCdt==1
                            ax = gca;
                            axPos = ax.Position;
                            axPos(1) = axPos(1)-.035;
                            axPos(3) = .03;
                            axes('Position',axPos);
                            hold on
                            
                            plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                            plot([0 0],[0 size(Profiles,1)],'k')
                            
                            tmp = abs(X_sort);
                            [~, idx] = min(tmp); %index of closest value
                            
                            MAX = ceil(max(abs(X_sort)));
                            
                            plot([MAX*-1 MAX],[idx idx],':k')
                            
                            axis([MAX*-1 MAX 1 size(Profiles,1)])
                            
                            t=ylabel(sprintf('%s\nPerc %s %s stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'Baseline', ...
                    ['Subj_' SubjID '_RasterBaseline_V_stim_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
    end
    
    
    
    %% Profiles stim = f(Percentile of A)
    close all
    Cdt =[1 1;1 2;1 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:3
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesCdtion{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesCdtion{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    All_Profiles_A{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_A{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                        
                        hold on
                        
                        colormap(ColorMap);
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-6 6])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if iROI==ROI_Groups{iRoiGrp}(1)
                            title([Name{Cdt(iCdt,2)} ' stimulus'])
                        end
                        
                        
                        if iCdt==1
                            ax = gca;
                            axPos = ax.Position;
                            axPos(1) = axPos(1)-.035;
                            axPos(3) = .03;
                            axes('Position',axPos);
                            hold on
                            
                            plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                            plot([0 0],[0 size(Profiles,1)],'k')
                            
                            tmp = abs(X_sort);
                            [~, idx] = min(tmp); %index of closest value
                            
                            MAX = ceil(max(abs(X_sort)));
                            
                            plot([MAX*-1 MAX],[idx idx],':k')
                            
                            axis([MAX*-1 MAX 1 size(Profiles,1)])
                            
                            t=ylabel(sprintf('%s\nPerc %s %s stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Percentile A stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'Baseline', ...
                    ['Subj_' SubjID '_RasterBaseline_A_stim_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
        
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles Att = f(Percentile of V)
    close all
    Cdt =[2 4; 2 1; 2 2; 2 3];
    Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
    NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            for iCdt = 1:4
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);

                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    if ROI(iROI).name(1)=='V'
                        Profiles=Profiles*-1;
                    end
                    
                    All_Profiles_AttfV{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_AttfV{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),4,iCdt+4*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-6 6])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        
                        if ROI(iROI).name(1)=='V'
                            t=title(NameSwitch{iCdt});
                            set(t,'fontsize',6)
                        else
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        if iCdt==1
                            ax = gca;
                            axPos = ax.Position;
                            axPos(1) = axPos(1)-.035;
                            axPos(3) = .03;
                            axes('Position',axPos);
                            hold on
                            
                            plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                            plot([0 0],[0 size(Profiles,1)],'k')
                            
                            tmp = abs(X_sort);
                            [~, idx] = min(tmp); %index of closest value
                            
                            MAX = ceil(max(abs(X_sort)));
                            
                            plot([MAX*-1 MAX],[idx idx],':k')
                            
                            axis([MAX*-1 MAX 1 size(Profiles,1)])
                            
                            t=ylabel(sprintf('%s\nPerc %s V stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Att = f(Percentile V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'Attention', ...
                    ['Subj_' SubjID '_RasterAttention_V_stim_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
        end
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles Att = f(Percentile of A)
    Cdt =[1 4; 1 1; 1 2; 1 3];
    Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
    NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:4
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end                    
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    if ROI(iROI).name(1)=='V'
                        Profiles=Profiles*-1;
                    end
                    
                    All_Profiles_AttfA{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_AttfA{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),4,iCdt+4*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-6 6])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        
                        if ROI(iROI).name(1)=='V'
                            t=title(NameSwitch{iCdt});
                            set(t,'fontsize',6)
                        else
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        if iCdt==1
                            ax = gca;
                            axPos = ax.Position;
                            axPos(1) = axPos(1)-.035;
                            axPos(3) = .03;
                            axes('Position',axPos);
                            hold on
                            
                            plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                            plot([0 0],[0 size(Profiles,1)],'k')
                            
                            tmp = abs(X_sort);
                            [~, idx] = min(tmp); %index of closest value
                            
                            MAX = ceil(max(abs(X_sort)));
                            
                            plot([MAX*-1 MAX],[idx idx],':k')
                            
                            axis([MAX*-1 MAX 1 size(Profiles,1)])
                            
                            t=ylabel(sprintf('%s\nPerc %s A stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Att = f(Percentile A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'Attention', ...
                    ['Subj_' SubjID '_RasterAttention_A_stim_' ToPlot{iToPlot}  '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles Att = f(Percentile of Att)
    close all
    Cdt =[4 4;1 1;1 1;3 3];
    Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
    NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:4
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    if ROI(iROI).name(1)=='V'
                        X=X*-1;
                    end
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    if ROI(iROI).name(1)=='V'
                        Profiles=Profiles*-1;
                    end
                    
                    All_Profiles_Att{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_Att{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),4,iCdt+4*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-4 4])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        
                        if ROI(iROI).name(1)=='V'
                            t=title(NameSwitch{iCdt});
                            set(t,'fontsize',6)
                        else
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        ax = gca;
                        axPos = ax.Position;
                        axPos(1) = axPos(1)-.032;
                        axPos(3) = .03;
                        axes('Position',axPos);
                        hold on
                        
                        plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                        plot([0 0],[0 size(Profiles,1)],'k')
                        
                        tmp = abs(X_sort);
                        [~, idx] = min(tmp); %index of closest value
                        
                        MAX = ceil(max(abs(X_sort)));
                        
                        plot([MAX*-1 MAX],[idx idx],':k')
                        
                        axis([MAX*-1 MAX 1 size(Profiles,1)])
                        
                        set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                            'ytick', linspace(1,size(Profiles,1),5),'yticklabel', [], ...
                            'ticklength', [0.01 0.01], 'fontsize', 6)
                        
                        
                        if iCdt==1
                            t=ylabel(sprintf('%s\nPerc Att contrast',...
                                strrep(ROI(iROI).name,'_','-')));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Att = f(Att) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'Attention', ...
                    ['Subj_' SubjID '_RasterAttention_Att_' ToPlot{iToPlot}  '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
        end
        
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles CrossMod = f(Percentile of V)
    close all
    Cdt =[2 1;2 2];
    Name = {'AV-A';'AV-V'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:2
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end                    
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    All_Profiles_Cross_V{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_Cross_V{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),2,iCdt+2*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-4 4])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if iROI == ROI_Groups{iRoiGrp}(1)
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        ax = gca;
                        axPos = ax.Position;
                        axPos(1) = axPos(1)-.032;
                        axPos(3) = .03;
                        axes('Position',axPos);
                        hold on
                        
                        plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                        plot([0 0],[0 size(Profiles,1)],'k')
                        
                        tmp = abs(X_sort);
                        [~, idx] = min(tmp); %index of closest value
                        
                        MAX = ceil(max(abs(X_sort)));
                        
                        plot([MAX*-1 MAX],[idx idx],':k')
                        
                        axis([MAX*-1 MAX 1 size(Profiles,1)])
                        
                        set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                            'ytick', linspace(1,size(Profiles,1),5),'yticklabel', [], ...
                            'ticklength', [0.01 0.01], 'fontsize', 10)
                        
                        
                        if iCdt==1
                            t=ylabel(sprintf('%s\nPerc %s V stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                        
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - CrossSensory = f(V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'CrossSensory', ...
                    ['Subj_' SubjID '_Raster_CrossSensory_V_stim_' ToPlot{iToPlot}   '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
        end
        
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles CrossMod = f(Percentile of A)
    close all
    Cdt =[1 1;1 2];
    Name = {'AV-A';'AV-V'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:2
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end                    
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    All_Profiles_Cross_A{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_Cross_A{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),2,iCdt+2*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-4 4])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if iROI == ROI_Groups{iRoiGrp}(1)
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        
                        ax = gca;
                        axPos = ax.Position;
                        axPos(1) = axPos(1)-.032;
                        axPos(3) = .03;
                        axes('Position',axPos);
                        hold on
                        
                        plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                        plot([0 0],[0 size(Profiles,1)],'k')
                        
                        tmp = abs(X_sort);
                        [~, idx] = min(tmp); %index of closest value
                        
                        MAX = ceil(max(abs(X_sort)));
                        
                        plot([MAX*-1 MAX],[idx idx],':k')
                        
                        axis([MAX*-1 MAX 1 size(Profiles,1)])
                        
                        set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                            'ytick', linspace(1,size(Profiles,1),5),'yticklabel', [], ...
                            'ticklength', [0.01 0.01], 'fontsize', 10)
                        
                        
                        if iCdt==1
                            t=ylabel(sprintf('%s\nPerc %s V stim',...
                                strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                            set(t,'fontsize',6)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - CrossSensory = f(A Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'CrossSensory', ...
                    ['Subj_' SubjID '_Raster_CrossSensory_A_stim_' ToPlot{iToPlot}  '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
    end
    
    clear Cdt Name ToPlot
    
    
    
    %% Profiles CrossMod = f(CrossMod)
    close all
    Cdt =[1 1;2 2];
    Name = {'AV-A';'AV-V'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:2
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end                    
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    All_Profiles_Cross{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_Cross{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),2,iCdt+2*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-4 4])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-3 3])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if iROI == ROI_Groups{iRoiGrp}(1)
                            t=title(Name{iCdt});
                            set(t,'fontsize',6)
                        end
                        
                        
                        ax = gca;
                        axPos = ax.Position;
                        axPos(1) = axPos(1)-.032;
                        axPos(3) = .03;
                        axes('Position',axPos);
                        hold on
                        
                        plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                        plot([0 0],[0 size(Profiles,1)],'k')
                        
                        tmp = abs(X_sort);
                        [~, idx] = min(tmp); %index of closest value
                        
                        MAX = ceil(max(abs(X_sort)));
                        
                        plot([MAX*-1 MAX],[idx idx],':k')
                        
                        axis([MAX*-1 MAX 1 size(Profiles,1)])
                        
                        set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                            'ytick', linspace(1,size(Profiles,1),5),'yticklabel', [], ...
                            'ticklength', [0.01 0.01], 'fontsize', 6)
                        
                        
                        if iCdt==1
                            t=ylabel(sprintf('%s\nPerc CrossSensory contrast',...
                                strrep(ROI(iROI).name,'_','-')));
                            set(t,'fontsize',10)
                            
                            set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                                'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                                'ticklength', [0.01 0.01], 'fontsize', 6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - CrossSensory = f(CrossSensory) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'CrossSensory', ...
                    ['Subj_' SubjID '_Raster_CrossSensory_' ToPlot{iToPlot}  '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
        
    end
    clear Cdt Name ToPlot
    
    
    
    %% Profiles Att = f(CrossMod)
    close all
    Cdt =[1 4;2 4];
    Name = {'AV-A';'AV-V'};
    NameAtt={'Att_A - Att_V'};
    NameSwitch={'Att_V - Att_A'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            if Print
                figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            end
            
            for iCdt = 1:2
                
                X_lh = nan(1,NbVertex(1));
                X_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,1));
                X_rh = nan(1,NbVertex(2));
                X_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,1));
                
                Profiles_lh = nan(NbVertex(1),6);
                Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1};
                
                Profiles_rh = nan(NbVertex(2),6);
                Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2};
                
                for iROI = ROI_Groups{iRoiGrp}
                    
                    X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                    [X_sort,I] = sort(X);
                    
                    Profiles = [...
                        Profiles_lh(ROI(iROI).VertOfInt{1},:) ; ...
                        Profiles_rh(ROI(iROI).VertOfInt{2},:)];
                    Profiles = Profiles(I,:);
                    
                    ToRemove = cat(3,isnan(Profiles), Profiles==0);
                    ToRemove = any(ToRemove,3);
                    ToRemove = any(ToRemove,2);
                    
                    Profiles(ToRemove,:)=[];
                    X_sort(ToRemove)=[];
                    
                    NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                    if numel(X)<NbBin
                        error('too many bins')
                    end                    
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    parfor iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort = X_sort_Perc;
                    Profiles = Profiles_Perc;
                    
                    if ROI(iROI).name(1)=='V'
                        Profiles=Profiles*-1;
                    end
                    
                    All_Profiles_Att_Cross{SubjInd,iToPlot,iCdt,iROI} = Profiles;
                    All_X_sort_Att_Cross{SubjInd,iToPlot,iCdt,iROI} = X_sort;
                    
                    if Print
                        %                     subplot(numel(ROI),2,iCdt+2*(iROI-1))
                        subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                        hold on
                        
                        colormap(ColorMap);
                        %                     imagesc(fliplr(Profiles), [-4 4])
                        imagesc(imgaussfilt(Profiles,[20 .001]), [-2 2])
                        
                        axis([0.5 6.5 0 size(Profiles,1)])
                        
                        DephLevels = round(linspace(100,0,8));
                        DephLevels([1;end]) = [];
                        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                            'ticklength', [0.01 0], 'fontsize', 4)
                        
                        t=xlabel('cortical depth');
                        set(t,'fontsize',6)
                        
                        if ROI(iROI).name(1)=='V'
                            t=title(NameSwitch);
                            set(t,'fontsize',6)
                        else
                            t=title(NameAtt);
                            set(t,'fontsize',6)
                        end
                        
                        
                        ax = gca;
                        axPos = ax.Position;
                        axPos(1) = axPos(1)-.032;
                        axPos(3) = .03;
                        axes('Position',axPos);
                        hold on
                        
                        plot(X_sort,1:numel(X_sort), 'b', 'linewidth',2)
                        plot([0 0],[0 size(Profiles,1)],'k')
                        
                        tmp = abs(X_sort);
                        [~, idx] = min(tmp); %index of closest value
                        
                        MAX = ceil(max(abs(X_sort)));
                        
                        plot([MAX*-1 MAX],[idx idx],':k')
                        
                        axis([MAX*-1 MAX 1 size(Profiles,1)])
                        
                        set(gca,'tickdir', 'out', 'xtick', [MAX*-1 MAX],'xticklabel', [MAX*-1 MAX], ...
                            'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
                            'ticklength', [0.01 0.01], 'fontsize', 6)
                        
                        
                        if iCdt==1
                            t=ylabel(sprintf('%s\nPerc %s',...
                                strrep(ROI(iROI).name,'_','-'), Name{Cdt(iCdt,1)}));
                            set(t,'fontsize',10)
                        else
                            t=ylabel(sprintf('Perc %s',...
                                Name{Cdt(iCdt,1)}));
                            set(t,'fontsize',6)
                        end
                    end
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
            if Print
                mtit(['Subject ' SubjID ' - Att = f(CrossSensory) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
                
                print(gcf, fullfile(FigureFolder,'CrossSensory', ...
                    ['Subj_' SubjID '_Raster_Att_CrossSensory_' ToPlot{iToPlot}   '_' ...
                    ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
                close all
            end
            
        end
    end
    
    clear Cdt Name ToPlot
    
    
    
    %%
    clear BetaCdtion BetaCrossSens BetaAtt FeaturesCdtion FeaturesAtt
    
end

cd(StartFolder)

save(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster.mat'), ...
    'ROI', ...
    'All_X_sort_V', 'All_X_sort_A', ...
    'All_X_sort_Att', 'All_X_sort_AttfV', 'All_X_sort_AttfA', ...
    'All_X_sort_Att_Cross',...
    'All_X_sort_Cross', 'All_X_sort_Cross_V', 'All_X_sort_Cross_A',...
    'All_Profiles_V', 'All_Profiles_A', ...
    'All_Profiles_Att', 'All_Profiles_AttfV', 'All_Profiles_AttfA', ...
    'All_Profiles_Att_Cross',...
    'All_Profiles_Cross', 'All_Profiles_Cross_V', 'All_Profiles_Cross_A')

return

%%
% Color map
X = 0:0.001:1;
R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
ColorMap = [R' G' B'];
clear R G B

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

StartFolder='/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

load(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster.mat'))

Subj2Include = true(11,1);
% Subj2Include([4 11]) = false;

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

ToPlot={'Cst','Lin'};

% ROI_Groups = {1:4;5:8;9:12;13:16;17:20;21:24;25:28};
% ROI_Groups_names = {'A1-PT-V123';'V123-ActDeact';'A1-PT-A-ActDeact';'V123-A-ActDeact';...
%     'A1-PT-V-ActDeact';'V123-V-ActDeact';'A1-PT-ActDeact'};
ROI_Groups = {1:4};
ROI_Groups_names = {'A1-PT-V123'};

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
% ToPermute = [];


%% Grp level  ; raster stim = f(Percentile of V)
close all
Cdt =[2 1;2 2;2 3];
Name={'A','V','AV'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_V,1)
                    X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-4 4])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI==ROI_Groups{iRoiGrp}(1)
                    title([Name{Cdt(iCdt,2)} ' stimulus'])
                end
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .15, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s %s stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_V_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')

    end
    
end


%% Grp level  ; raster stim = f(Percentile of A)
Cdt =[1 1;1 2;1 3];
Name={'A','V','AV'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_A,1)
                    X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-4 4])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .15, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s %s stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Percentile A stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_A_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster Att = f(Percentile of V)
close all
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:4
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_A,1)
                    X_sort(iSubj,:) = All_X_sort_AttfV{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                
                subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch{iCdt});
                    set(t,'fontsize',6)
                else
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .1, .1, .04, .04, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf(sprintf('%s\nPerc %s V stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['GrpLvl_Raster_Attention_V_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
    end
end


%% Grp level  ; raster Att = f(Percentile of A)
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:4
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_A,1)
                    X_sort(iSubj,:) = All_X_sort_AttfA{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch{iCdt});
                    set(t,'fontsize',6)
                else
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .1, .1, .04, .04, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s A stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['GrpLvl_Raster_Attention_A_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
        
        
    end
end



%% Grp level  ; raster CrossMod = f(Percentile of V stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);

                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI == ROI_Groups{iRoiGrp}(1)
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);

                YLabel = sprintf('%s\nPerc %s V stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                
                clear X
                
            end
            
        end
        
        mtit(['CrossMod = f(Percentile V Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_CrossMod_Stim_V_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster CrossMod = f(Percentile of A stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);                
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI == ROI_Groups{iRoiGrp}(1)
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                YLabel = sprintf('%s\nPerc %s A stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                
                clear X
                
            end
            
        end
        
        mtit(['CrossMod = f(Percentile A Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_CrossMod_Stim_A_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
        
        
    end
end



%% Grp level  ; raster Att = f(Percentile of CrossMod)
close all
Cdt =[1 4;2 4];
Name = {'AV-A';'AV-V'};
NameAtt={'Att_A - Att_V'};
NameSwitch={'Att_V - Att_A'};
ToPlot={'Cst','Lin'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    X_sort(iSubj,:) = All_X_sort_Att_Cross{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Att_Cross{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);                
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-1 1])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch);
                    set(t,'fontsize',6)
                else
                    t=title(NameAtt);
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s',...
                        strrep(ROI(iROI).name,'_','-'), Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                else
                    YLabel = sprintf('Perc %s',...
                        Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                    
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile CrossMod) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_Att_CrossMod_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '.tif']), '-dtiff')
        
        
    end
    
end

