clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces'); %#ok<*NASGU>

SubjectList = [...
%     '02';...
%     '03';...
    '04';...
%     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
%     '13';...
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

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')

%%
for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder, 'BetaMapping','8Surf');
    
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
        
        
        
        FeatureSaveFile = fullfile(Data_Folder,[ 'Subj_' SubjID '_features_' HsSufix 'hs_' ...
            num2str(NbLayers) '_surf.mat']);
        
%         InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
%             ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
%         
%         [Vertex,~,~] = read_vtk(InfSurfFile, 0, 1);
%         NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData; %#ok<*SAGROW>
        else
            Betas = dir(fullfile(Data_Folder, ['Beta*' HsSufix 'cr.vtk']));
            
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
        for CondInd = 1:length(Conditions_Names)% For each Condition
            
            Beta2Sel = [];
            Beta2Sel2 = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd} ' - Block ' num2str(BlockInd) '*bf(1)']))];
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            
            Features = AllMapping(:,2:end-1,Beta2Sel);
            FeaturesCdtion{CondInd,hs,1} = mean(Features(:,:,1:(size(Features,3)/2)),3);
            FeaturesCdtion{CondInd,hs,2} = mean(Features(:,:,(1+(size(Features,3)/2)):end),3);
            
            
            clear Features Beta2Sel Beta2Sel2 X Y Mapping
            
        end
        
        
        %% Run GLMs for A, V and AV for both attention
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
            
            X=repmat(DesMat,size(Features,3)/2,1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            Y1 = Y(1:size(Y,1)/2,:);
            Y2 = Y((1+size(Y,1)/2):end,:);
            BetaCdtion{hs}(:,:,CondInd,1) = pinv(X)*Y1;
            BetaCdtion{hs}(:,:,CondInd,2) = pinv(X)*Y2;
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel Beta2Sel2 X Y Mapping
            
        end
        
        clear CondInd AttCondInd CrossSensCondInd CrossSensCondNames AttCondNames BlockInd CondInd
        
    end
    
    
%     if any(NbVertex ~= NbVertices)
%         %         NbVertex
%         %         NbVertices
%         error('The number of vertices does not match.')
%     end
    
    close all
    
    
    
    %% Profiles stim = f(Percentile of V)
    close all
    Cdt =[2 1;2 2;2 3;2 4;2 5;2 6];
    Name={'A','V','AV','A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        close all
        
        for iROI = 1:4 %numel(ROI)
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            fprintf('    %s\n',ROI(iROI).name)
            
            for iCdt = 1:6
                
                CV_cdt = [2 1];
                
                clear X_sort_CV Profiles_CV
                
                for iCV =1:2
                    X_lh = nan(1,NbVertex(1));
                    X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                    X_rh = nan(1,NbVertex(2));
                    X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                    
                    Profiles_lh = nan(NbVertex(1),6);
                    Profiles_lh(VertexWithDataHS{1},:) = FeaturesCdtion{Cdt(iCdt,2),1,CV_cdt(iCV)};
                    
                    Profiles_rh = nan(NbVertex(2),6);
                    Profiles_rh(VertexWithDataHS{2},:) = FeaturesCdtion{Cdt(iCdt,2),2,CV_cdt(iCV)};
                    
                    
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
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    for iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    clear X_sort Profiles
                    
                    X_sort_CV(:,iCV) = X_sort_Perc;
                    Profiles_CV(:,:,iCV) = Profiles_Perc;
                end
                
                All_Profiles_V{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                All_X_sort_V{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);
                
                clear X
                
            end
            
            clear X_lh Y_rh
            
        end
        
        
    end
    
    
    
    %% Profiles stim = f(Percentile of A)
    Cdt =[1 1;1 2;1 3;1 4;1 5;1 6];
    Name={'A','V','AV','A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        close all
        
        for iROI = 1:4 %numel(ROI)
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            fprintf('    %s\n',ROI(iROI).name)
            
            for iCdt = 1:6
                
                CV_cdt = [2 1];
                
                clear X_sort_CV Profiles_CV
                
                for iCV =1:2
                    X_lh = nan(1,NbVertex(1));
                    X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                    X_rh = nan(1,NbVertex(2));
                    X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                    
                    Profiles_lh = nan(NbVertex(1),6);
                    Profiles_lh(VertexWithDataHS{1},:) = FeaturesCdtion{Cdt(iCdt,2),1,CV_cdt(iCV)};
                    
                    Profiles_rh = nan(NbVertex(2),6);
                    Profiles_rh(VertexWithDataHS{2},:) = FeaturesCdtion{Cdt(iCdt,2),2,CV_cdt(iCV)};
                    
                    
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
                    
                    IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                    
                    clear X_sort_Perc Profiles_Perc
                    
                    for iPerc = 2:numel(IdxToAvg)
                        X_sort_Perc(iPerc-1) = mean(X_sort(IdxToAvg((iPerc-1):iPerc)));
                        Profiles_Perc(iPerc-1,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:));
                    end
                    
                    X_sort_CV(:,iCV) = X_sort_Perc;
                    Profiles_CV(:,:,iCV) = Profiles_Perc;
                end
                
                All_Profiles_A{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                All_X_sort_A{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);
                
                clear X
                
            end
            
            clear X_lh Y_rh
            
        end
        
        
    end
    
    
    %%
    clear BetaCdtion BetaCrossSens BetaAtt FeaturesCdtion FeaturesAtt
    
end

cd(StartFolder)

save(fullfile(StartFolder,'Results','Profiles','Surfaces','RasterAllCdt_Lin_CV.mat'), ...
    'ROI', ...
    'All_X_sort_V', 'All_X_sort_A', ...
    'All_Profiles_V', 'All_Profiles_A')

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

load(fullfile(StartFolder,'Results','Profiles','Surfaces','RasterAllCdt_CV.mat'))

Subj2Include = true(11,1);
% Subj2Include([4 11]) = false;

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

ToPlot={'Cst','Lin'};

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')


CLIM = [-5 5];


%% Grp level  ; raster stim = f(Percentile of V)
close all
Cdt =[2 1;2 2;2 3;2 4;2 5;2 6];
Name={'A','V','AV','A','V','AV'};
ToPlot={'Cst','Lin'};

for iToPlot = 1%:numel(ToPlot)
    
    close all
    
    for iROI = 1:4 %:numel(ROI)
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        fprintf('    %s\n',ROI(iROI).name)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
        
        for iCdt = 1:6
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_V,1)
                Sorting_Raster(:,:,iSubj) = (All_Profiles_V{iSubj,iToPlot,2,iROI}+All_Profiles_V{iSubj,iToPlot,5,iROI})/2;
                X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            MeanProfiles = mean(Profiles,3);
            
            
            subplot(2,3,iCdt)
            
            colormap(ColorMap);
            imagesc(flipud(imgaussfilt(MeanProfiles,[10 .001])), CLIM)
            %             imagesc(MeanProfiles, [-5 5])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t = title(strrep(Conditions_Names{iCdt},'ention', ''));
            set(t,'fontsize',10)
            
            
            ax = gca;
            axes('Position',ax.Position);
            hold on
            errorbar(1:6,mean(MeanProfiles),nansem(squeeze(mean(Profiles)),2), 'k', 'linewidth', 2)
            for iSubj=1:size(Profiles,3)
                plot(1:6,mean(Profiles(:,:,iSubj)), ':k', 'linewidth', .5)
            end
            plot([1 6],[0 0], '--k')
            
            axis([0.5 6.5 -3 6])
            
            DephLevels = round(linspace(100,0,8));
            DephLevels([1;end]) = [];
            set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                'YAxisLocation','right', 'ytick', -10:10,'yticklabel', -10:10, ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            
            if iCdt==1 ||  iCdt==4
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s %s stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
            end
            
            if iCdt==3 ||  iCdt==6
                PlotColorBar(ax, ColorMap, CLIM)
            end
        end
        
        mtit([strrep(ROI(iROI).name,'_',' ') ' - Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_RasterAllCdt_Baseline_V_stim_' ToPlot{iToPlot} '_' ROI(iROI).name '_CV.tiff']), '-dtiff')
        
    end
    
end


%% Grp level  ; raster stim = f(Percentile of A)
close all
Cdt =[1 1;1 2;1 3;1 4;1 5;1 6];
Name={'A','V','AV','A','V','AV'};
ToPlot={'Cst','Lin'};

for iToPlot = 1%:numel(ToPlot)
    
    close all
    
    for iROI = 1:4 %numel(ROI)
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        fprintf('    %s\n',ROI(iROI).name)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
        
        for iCdt = 1:6
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_V,1)
                Sorting_Raster(:,:,iSubj) = (All_Profiles_A{iSubj,iToPlot,1,iROI}+All_Profiles_V{iSubj,iToPlot,4,iROI})/2;
                X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
            end
            
            Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            MeanProfiles = mean(Profiles,3);
            
            
            subplot(2,3,iCdt)
            
            colormap(ColorMap);
            imagesc(flipud(imgaussfilt(MeanProfiles,[10 .001])), CLIM)
            %             imagesc(MeanProfiles, [-5 5])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t = title(strrep(Conditions_Names{iCdt},'ention', ''));
            set(t,'fontsize',10)
            
            
            ax = gca;
            axes('Position',ax.Position);
            hold on
            errorbar(1:6,mean(MeanProfiles),nansem(squeeze(mean(Profiles)),2), 'k', 'linewidth', 2)
            for iSubj=1:size(Profiles,3)
                plot(1:6,mean(Profiles(:,:,iSubj)), ':k', 'linewidth', .5)
            end
            plot([1 6],[0 0], '--k')
            
            axis([0.5 6.5 -3 6])
            
            DephLevels = round(linspace(100,0,8));
            DephLevels([1;end]) = [];
            set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                'YAxisLocation','right', 'ytick', -10:10,'yticklabel', -10:10, ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            
            if iCdt==1 ||  iCdt==4
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s %s stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
            end
            
            if iCdt==3 ||  iCdt==6
                PlotColorBar(ax, ColorMap, CLIM)
            end
            
        end
        
        mtit([strrep(ROI(iROI).name,'_',' ') ' - Percentile A stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_RasterAllCdt_Baseline_A_stim_' ToPlot{iToPlot} '_' ROI(iROI).name '_CV.tiff']), '-dtiff')
        
    end
    
end
