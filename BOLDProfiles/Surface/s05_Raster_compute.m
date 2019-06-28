% computes all rasters


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
        
        % Load data or extract them
        [AllMapping, Face, Vertex, VertexWithData] = ...
            load_vtk_beta_maps(FeatureSaveFile, Data_Folder, HsSufix, NbVertices, NbLayers);
        
        VertexWithDataHS{hs} = VertexWithData;
        
        
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
            FeaturesCdtion{CondInd,hs,1} = mean(Features(:,:,1:(size(Features,3)/2)),3);
            FeaturesCdtion{CondInd,hs,2} = mean(Features(:,:,(1+(size(Features,3)/2)):end),3);
            
            X=repmat(DesMat,size(Features,3)/2,1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            Y1 = Y(1:size(Y,1)/2,:);
            Y2 = Y((1+size(Y,1)/2):end,:);
            BetaCdtion{hs}(:,:,CondInd,1) = pinv(X)*Y1;
            BetaCdtion{hs}(:,:,CondInd,2) = pinv(X)*Y2;
            
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
            
            FeaturesAtt{AttCondInd,hs,1} = mean(Features(:,:,1:(size(Features,3)/2)),3);
            FeaturesAtt{AttCondInd,hs,2} = mean(Features(:,:,(1+(size(Features,3)/2)):end),3);
            
            FeaturesAll(:,:,:,AttCondInd) = Features;
            
            X=repmat(DesMat,size(Features,3)/2,1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            Y1 = Y(1:size(Y,1)/2,:);
            Y2 = Y((1+size(Y,1)/2):end,:);
            BetaAtt{hs}(:,:,AttCondInd,1) = pinv(X)*Y1;
            BetaAtt{hs}(:,:,AttCondInd,2) = pinv(X)*Y2;
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        Features = mean(FeaturesAll,4);
        FeaturesAtt{size(Cond2Contrast,1)+1,hs,1} = mean(Features(:,:,1:(size(Features,3)/2)),3);
        FeaturesAtt{size(Cond2Contrast,1)+1,hs,2} = mean(Features(:,:,(1+(size(Features,3)/2)):end),3);
        
        X=repmat(DesMat,size(Features,3)/2,1);
        
        Y = shiftdim(Features,1);
        Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
        
        Y1 = Y(1:size(Y,1)/2,:);
        Y2 = Y((1+size(Y,1)/2):end,:);
        BetaAtt{hs}(:,:,end+1,1) = pinv(X)*Y1;
        BetaAtt{hs}(:,:,end,2) = pinv(X)*Y2;
        
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
            
            FeaturesCrossMod{CrossSensCondInd,hs,1} = mean(Features(:,:,1:(size(Features,3)/2)),3);
            FeaturesCrossMod{CrossSensCondInd,hs,2} = mean(Features(:,:,(1+(size(Features,3)/2)):end),3);
            
            X=repmat(DesMat,size(Features,3)/2,1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            Y1 = Y(1:size(Y,1)/2,:);
            Y2 = Y((1+size(Y,1)/2):end,:);
            BetaCrossSens{hs}(:,:,CrossSensCondInd,1) = pinv(X)*Y1;
            BetaCrossSens{hs}(:,:,CrossSensCondInd,2) = pinv(X)*Y2;
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        clear CondInd AttCondInd CrossSensCondInd CrossSensCondNames AttCondNames BlockInd CondInd
        
    end
    

    %% Profiles stim = f(Percentile of V)
    Cdt =[2 1;2 2;2 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:3
                    
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
                        
                        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                        
                        IdxToAvg = floor(linspace(1,numel(X_sort),NbBin+1));
                        
                        clear X_sort_Perc Profiles_Perc
                        
                        parfor iPerc = 2:numel(IdxToAvg)
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
    end
    
    
    
    %% Profiles stim = f(Percentile of A)
    close all
    Cdt =[1 1;1 2;1 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iRoiGrp = 1:numel(ROI_Groups)
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:3
                    
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:4
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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
                        
                        if ROI(iROI).name(1)=='V'
                            Profiles_Perc=Profiles_Perc*-1;
                        end

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_AttfV{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_AttfV{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);   
                                      
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:4
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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
                        
                        if ROI(iROI).name(1)=='V'
                            Profiles_Perc=Profiles_Perc*-1;
                        end

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_AttfA{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_AttfA{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);   
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:4
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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
                        
                        if ROI(iROI).name(1)=='V'
                            Profiles_Perc=Profiles_Perc*-1;
                        end

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_Att{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_Att{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);                    
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:2
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_Cross_V{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_Cross_V{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);                        
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:2
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_Cross_A{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_Cross_A{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);                      
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:2
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesCrossMod{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesCrossMod{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_Cross{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_Cross{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);
                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
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
            
            for iROI = ROI_Groups{iRoiGrp}
                
                for iCdt = 1:2
                    
                    CV_cdt = [2 1];
                    
                    clear X_sort_CV Profiles_CV
                    
                    for iCV =1:2
                        
                        X_lh = nan(1,NbVertex(1));
                        X_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,1),iCV);
                        X_rh = nan(1,NbVertex(2));
                        X_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,1),iCV);
                        
                        Profiles_lh = nan(NbVertex(1),6);
                        Profiles_lh(VertexWithDataHS{1},:) = FeaturesAtt{Cdt(iCdt,2),1,CV_cdt(iCV)};
                        
                        Profiles_rh = nan(NbVertex(2),6);
                        Profiles_rh(VertexWithDataHS{2},:) = FeaturesAtt{Cdt(iCdt,2),2,CV_cdt(iCV)};
                        
                        
                        
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
                        
                        if ROI(iROI).name(1)=='V'
                            Profiles_Perc=Profiles_Perc*-1;
                        end

                        X_sort_CV(:,iCV) = X_sort_Perc;
                        Profiles_CV(:,:,iCV) = Profiles_Perc;
                        
                    end
                    
                    All_Profiles_Att_Cross{SubjInd,iToPlot,iCdt,iROI} = mean(Profiles_CV,3);
                    All_X_sort_Att_Cross{SubjInd,iToPlot,iCdt,iROI} = mean(X_sort_CV,2);

                    
                    clear X
                    
                end
                
                clear X_lh Y_rh
                
            end
            
        end
    end
    
    clear Cdt Name ToPlot
    
    
    
    %%
    clear BetaCdtion BetaCrossSens BetaAtt FeaturesCdtion FeaturesAtt
    
end

cd(StartFolder)

save(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster_CV.mat'), ...
    'ROI', ...
    'All_X_sort_V', 'All_X_sort_A', ...
    'All_X_sort_Att', 'All_X_sort_AttfV', 'All_X_sort_AttfA', ...
    'All_X_sort_Att_Cross',...
    'All_X_sort_Cross', 'All_X_sort_Cross_V', 'All_X_sort_Cross_A',...
    'All_Profiles_V', 'All_Profiles_A', ...
    'All_Profiles_Att', 'All_Profiles_AttfV', 'All_Profiles_AttfA', ...
    'All_Profiles_Att_Cross',...
    'All_Profiles_Cross', 'All_Profiles_Cross_V', 'All_Profiles_Cross_A')


