clear; close all; clc;

NbLayers = 6;

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ExtFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/';

ROIs = {...
    'TE', {'TE1.0','TE1.1','TE1.2'}
    'pSTG', {'pSTG'}; ...
    'V1', {'V1'}; ...
    'V2-3', {'V2', 'V2'}; ...
    };

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));
[~,~,~] = mkdir(FigureFolder);

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

NbSubject = size(SubjectList,1);

for iROI = 1:length(ROIs)
    AllSubjects_Data(iROI) = struct(...
        'name', ROIs{iROI}, ...
        'DATA', nan(NbLayers,6,NbSubject), ...
        'VoxelCount', nan(NbSubject,NbLayers+1),...
        'ROI_Coverage', nan(NbSubject,2),...
        'Differential', struct('Blocks', struct('DATA', cell(1))));
end


%% Gets data for each subject
for SubjInd = 1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    SubjectFolder = fullfile(ExtFolder, ['Subject_' SubjID]);
    AnalysisFolder = fullfile(SubjectFolder, 'BetaMapping', [num2str(NbLayers+2) 'Surf']);
    
    
    for iROI=1:size(ROIs,1)
        
        %% Get the data
        
        if exist(fullfile(SubjectFolder, 'BetaMapping', [num2str(NbLayers+2) 'Surf'], ...
                strcat('Data_', ROIs{iROI,1}, '_', num2str(NbLayers+2), '_Surf.mat')), 'file')
            load(fullfile(SubjectFolder, 'BetaMapping', [num2str(NbLayers+2) 'Surf'], ...
                strcat('Data_', ROIs{iROI,1}, '_', num2str(NbLayers+2), '_Surf.mat')))
        else
            
            VertCount = [0 0];
            DATA=[];
            
            for iSubROI=1:numel(ROIs{iROI,2})
                for i=1:2
                    if i==1
                        File2Load = fullfile(AnalysisFolder, strcat('Features_', ...
                            ROIs{iROI,2}{iSubROI}, '_lhs_', num2str(NbLayers+2), '_Surf.mat'));
                    else
                        File2Load = fullfile(AnalysisFolder, strcat('Features_', ...
                            ROIs{iROI,2}{iSubROI}, '_rhs_', num2str(NbLayers+2), '_Surf.mat'));
                    end
                    
                    load(File2Load)
                    
                    VertCount(i) = VertCount(i)+size(Features,2);
                    DATA = [DATA Features];
                    
                    clear Features FilesCell File2Load
                end
                clear i
            end
            clear iSubROI
            
            AllSubjects_Data(iROI).VerticesCount(SubjInd,1:2) = VertCount; clear VertCount
            Features = DATA; clear DATA
            
            LogFeat = ~any(isnan(Features));
            FeaturesLayers = repmat((2+NbLayers):-1:1, 1, size(Features,2)/(2+NbLayers));
            
            
            %% Get which line corresponds to which condition
            load(fullfile(SubjectFolder, 'FFX_Block', 'SPM.mat'))
            
            BetaNames = char(SPM.xX.name');
            
            BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
                strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
                strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
                strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
                strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
                2);
            
            BetaOfInterest = find(BetaOfInterest);
            
            BetaNames(:,1:6)=[];
            
            for iCond=1:numel(Conditions_Names)
                
                tmp=BetaNames(BetaOfInterest,1:length(Conditions_Names{iCond}));
                CondLines(:,iCond) = find(strcmp(Conditions_Names{iCond}, cellstr(tmp)));
                
                clear tmp
            end
            
            clear iCond SPM BetaOfInterest BetaNames
            
            
            
            %% Averages across blocks and voxels
            for iCond = 1:numel(Conditions_Names) % For each Condition
                
                Beta2Sel = CondLines(:,iCond);
                
                for BlockInd = 1:size(Beta2Sel,1) % For each Block
                    
                    Ind = NbLayers:-1:1;
                    
                    for LayerInd = 1:NbLayers % Averages over voxels of a given layer
                        
                        Data_ROI.LayerMean(Ind(LayerInd),BlockInd,iCond) = ...
                            nanmean(Features(Beta2Sel(BlockInd),FeaturesLayers==LayerInd+1));
                        Data_ROI.LayerSTD(Ind(LayerInd),BlockInd,iCond) = ...
                            nanstd(Features(Beta2Sel(BlockInd),FeaturesLayers==LayerInd+1));
                        Data_ROI.LayerSEM(Ind(LayerInd),BlockInd,iCond)  = ...
                            nansem(Features(Beta2Sel(BlockInd),FeaturesLayers==LayerInd+1));
                    end
                    clear LayerInd
                    
                end
                clear BlockInd
                
            end
            
            Data_ROI.MEAN=squeeze(nanmean(Data_ROI.LayerMean,2));
            Data_ROI.STD=squeeze(nanstd(Data_ROI.LayerMean,2));
            Data_ROI.SEM=squeeze(nansem(Data_ROI.LayerMean,2));
            
            save(fullfile(SubjectFolder, 'BetaMapping', [num2str(NbLayers+2) 'Surf'], ...
                strcat('Data_', ROIs{iROI,1}, '_', num2str(NbLayers+2), '_Surf.mat')), 'Data_ROI')
        
        end       
        
        
        %% Compute main effect for each block
        A=Data_ROI.LayerMean(1:end,:,:);
        
        Blocks = A(:,:,1:3);
        Blocks(:,:,end+1) = nanmean(A(:,:,1:3),3);
        Blocks(:,:,end+1:end+3) = A(:,:,4:6);
        Blocks(:,:,end+1) = nanmean(A(:,:,4:6),3);
        Blocks(:,:,end+1) = nanmean(A(:,:,[1 4]),3);
        Blocks(:,:,end+1) = nanmean(A(:,:,[2 5]),3);
        Blocks(:,:,end+1) = nanmean(A(:,:,[3 6]),3);
        Blocks(:,:,end+1) = nan(size(Blocks(:,:,end)));
        
        AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = Blocks;
        
        clear A Blocks
        
        %% Compute differential for each block
        A=Data_ROI.LayerMean(1:end,:,:);
        
        if any(A(:)==0)
            A(A==0)=NaN;
        end
        AllSubjects_Data(iROI).DATA(:,:,SubjInd) = squeeze(nanmean(A,2));
        
        A(:,:,end+1:end+3) = A(:,:,4:6)-A(:,:,1:3);
        
        
        for j=1:3
            TEMP(:,:,j) = A(:,:,j*3) - ( A(:,:,j*3-1)+A(:,:,j*3-2) ); %#ok<*SAGROW>
        end
        
        Blocks = A(:,:,1:3);
        Blocks(:,:,end+1) = TEMP(:,:,1);
        Blocks(:,:,end+1:end+3) = A(:,:,4:6);
        Blocks(:,:,end+1) = TEMP(:,:,2);
        Blocks(:,:,end+1:end+3) = A(:,:,7:9);
        Blocks(:,:,end+1) = TEMP(:,:,3);
        
        MainEffects{1} = mean(cat(3,Blocks(:,:,4),Blocks(:,:,8)),3);
        MainEffects{2} = mean(cat(3,Blocks(:,:,9),Blocks(:,:,10),Blocks(:,:,11)),3);
        
        AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjInd} = Blocks;
        AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
        
        
        clear A TEMP Blocks temp
        
        %% Compute AV-V and AV-A for each block
        A=Data_ROI.LayerMean(1:end,:,:);
        
        if any(A(:)==0)
            A(A==0)=NaN;
        end
        AllSubjects_Data(iROI).DATA(:,:,SubjInd) = squeeze(nanmean(A,2));
        
        Blocks = A;
        Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,1), A(:,:,6)-A(:,:,4)),3);
        Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,2), A(:,:,6)-A(:,:,5)),3);
        
        AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjInd} = Blocks;
        
        clear Data_ROIs SubROI_Ind Data_ROI A TEMP Blocks temp
        
    end
end

clear SubjInd ROI_Ind

cd(StartFolder)


%% Averages over subjects
for iROI=1:length(AllSubjects_Data)
    
    Include = true(NbSubject,1);
    
    AllSubjects_Data(iROI).Include=Include;
    
    AllSubjects_Data(iROI).MEAN = nanmean(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).STD = nanstd(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).SEM = nansem(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    
    % Compute differential
    A = AllSubjects_Data(iROI).DATA(:,:,Include);
    A = [A A(:,4:6,:)-A(:,1:3,:)]; %#ok<AGROW>
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    AllSubjects_Data(iROI).Differential.DATA = A;
    AllSubjects_Data(iROI).Differential.MEAN = nanmean(A,3);
    AllSubjects_Data(iROI).Differential.STD = nanstd(A,3);
    AllSubjects_Data(iROI).Differential.SEM = nansem(A,3);
    
    B(:,1,:) = nanmean(A(:,[4 8],:),2);
    B(:,2,:) = nanmean(A(:,9:11,:),2);
    
    AllSubjects_Data(iROI).Differential.MainEffect.DATA = B;
    AllSubjects_Data(iROI).Differential.MainEffect.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).Differential.MainEffect.STD = nanstd(B,3);
    AllSubjects_Data(iROI).Differential.MainEffect.SEM = nansem(B,3);
    
    clear A B
    
    % Compute main effects
    A = AllSubjects_Data(iROI).DATA(:,:,Include);
    TEMP(:,1:3,:) = A(:,1:3,:);
    TEMP(:,end+1,:) = nanmean(A(:,1:3,:),2);
    TEMP(:,end+1:end+3,:) = A(:,4:6,:);
    TEMP(:,end+1,:) = nanmean(A(:,4:6,:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[1 4],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[2 5],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[3 6],:),2);
    TEMP(:,end+1,:) = nan(size(TEMP(:,end,:)));
    
    AllSubjects_Data(iROI).MainEffects.DATA = TEMP;
    AllSubjects_Data(iROI).MainEffects.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).MainEffects.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).MainEffects.SEM = nansem(TEMP,3);
    
    clear A TEMP
    
    % Compute Bi - Uni
    A = AllSubjects_Data(iROI).DATA(:,:,Include);
    TEMP = A(:,1:6,:);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
    
    AllSubjects_Data(iROI).BiVSUni.DATA = TEMP;
    AllSubjects_Data(iROI).BiVSUni.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).BiVSUni.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).BiVSUni.SEM = nansem(TEMP,3);
    
    clear A TEMP
    
end




%% Computes betas
DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [DesMat' (DesMat.^2)' ones(NbLayers,1)];
DesMat = [DesMat' ones(NbLayers,1)];

for iROI=1:length(AllSubjects_Data)
    
    Name=AllSubjects_Data(iROI).name;
    
    fprintf('\nComputing betas for ROI %s\n', Name)
    
    SubjToInclude = find(AllSubjects_Data(iROI).Include);
    
    AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA = ...
        nan(size(DesMat,2),2,size(SubjectList,1));
    
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    
    for i=1:3
        
        for SubjInd = 1:length(SubjToInclude)
            
            switch i
                case 1
                    Blocks = AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjToInclude(SubjInd)};
                case 2
                    Blocks = AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjToInclude(SubjInd)};
                case 3
                    Blocks = AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjToInclude(SubjInd)};
            end
            
            if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
                
                for CondInd = 1:size(Blocks,3);
                    
                    Y = flipud(Blocks(:,:,CondInd));
                    if any(isnan(Y(:)))
                        [~,y]=find(isnan(Y));
                        y=unique(y);
                        Y(:,y)=[];
                        clear y
                    end
                    
                    if isempty(Y)
                        B=nan(1,size(DesMat,2));
                    else
                        X=repmat(DesMat,size(Y,2),1);
                        Y=Y(:);
                        [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
                    end
                    
                    switch i
                        case 1
                            AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 2
                            AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 3
                            AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                    end
                    clear X Y B
                    
                end
                
                
            end
            
            % Betas for main effect of differential
            if i==1
                for j=1:numel(AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)})
                    Blocks = AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)}{j};
                    if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
                        Y = flipud(Blocks(:,:));
                        if any(isnan(Y(:)))
                            [~,y]=find(isnan(Y));
                            y=unique(y);
                            Y(:,y)=[];
                            clear y
                        end
                        
                        if isempty(Y)
                            B=nan(1,size(DesMat,2));
                        else
                            X=repmat(DesMat,size(Y,2),1);
                            Y=Y(:);
                            [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
                        end
                        
                        AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA(:,j,SubjToInclude(SubjInd))=B;
                        
                        clear Y B
                    end
                    clear X
                end
            end
            
            
        end
        
    end
    
    % Main Effects average betas
    tmp=AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).MainEffects.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Differentials average betas
    tmp=AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).Differential.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).Differential.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).Differential.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).Differential.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Differentials main effects average betas
    tmp=AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA;
    
    AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Bi VS Uni average betas
    tmp=AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).BiVSUni.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    clear tmp H P
    
end


%% Saves
fprintf('\nSaving\n')

cd(FigureFolder)
save(strcat('Data_Block_Surf_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'))

cd(StartFolder)