clear; close all; clc;

NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

InclusionThreshold = .5;

ANTs = 0;
Median = 1;

StartFolder=fullfile(pwd, '..','..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ROIs = {...
    'A1_surf';...
    'PT_surf_thres';...
    'V1_surf_thres';...
    'V2-3_surf_thres';...

% 'V1_A_Deact';...
% 'V1_V_AV_Deact';...
% 'V1_V_AV_Act';...
% 'V2-3_A_Deact';...
% 'V2-3_V_AV_Deact';...
% 'V2-3_V_AV_Act';...
    };


FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles', 'Volumes', strcat(num2str(NbLayers), '_layers'));
[~,~,~] = mkdir(FigureFolder);

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

NbSubject = size(SubjectList,1);

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


%%
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
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    AnalysisFolder = fullfile(SubjectFolder,'Results', 'Profiles', 'Volumes', TargetLayerFile);
    
    for iROI=1:size(ROIs,1)
        
        if ANTs
            File2Load = fullfile(AnalysisFolder, strcat('Data_Block_', ...
                AllSubjects_Data(iROI).name, '_', TargetLayerFile, '_ANTs.mat'));
        else
            File2Load = fullfile(AnalysisFolder, strcat('Data_Block_', ...
                AllSubjects_Data(iROI).name, '_', TargetLayerFile, '.mat')); %#ok<*UNRCH>
        end
        
        clear Data_ROI
        
        if exist(File2Load, 'file')
            
            load(File2Load, 'Data_ROI')
            
            AllSubjects_Data(iROI).ROI_Coverage(SubjInd,:) = Data_ROI.ROI_Coverage;
            
            if any(Data_ROI.Voxelcount<50)
                
                AllSubjects_Data(iROI).VoxelCount(SubjInd,1:NbLayers+1) = nan(1,NbLayers+1);
                
                AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = nan(NbLayers,12,2);
                AllSubjects_Data(iROI).Interaction.Blocks.DATA{SubjInd} = nan(NbLayers,12,3);
                AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjInd} = nan(NbLayers,12,2);
                AllSubjects_Data(iROI).InteractionRestrict.Blocks.DATA{SubjInd} = nan(NbLayers,12,3);
%                 AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjInd} = nan(NbLayers,12,2);
%                 AllSubjects_Data(iROI).InteractionSkim.Blocks.DATA{SubjInd} = nan(NbLayers,12,3);
%                 AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjInd} = nan(NbLayers,12,2);
%                 AllSubjects_Data(iROI).InteractionPSC.Blocks.DATA{SubjInd} = nan(NbLayers,12,3);
                
            else
                
                AllSubjects_Data(iROI).VoxelCount(SubjInd,1:NbLayers+1) = Data_ROI.Voxelcount;
                
                
                
                for i=1:2
                    if Median
                        switch i
                            case 1
                                AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEDIAN(2:NbLayers+1,:);
                                A=Data_ROI.LayerMedian(2:NbLayers+1,:,:);
                            case 2
                                AllSubjects_Data(iROI).DATARestrict(:,:,SubjInd) = Data_ROI.MEDIANRestrict(2:NbLayers+1,:);
                                A=Data_ROI.LayerMedianRestrict(2:NbLayers+1,:,:);
                            case 3
%                                 AllSubjects_Data(iROI).DATASkim(:,:,SubjInd) = Data_ROI.MEDIANSkim(2:NbLayers+1,:);
%                                 A=Data_ROI.LayerMedianSkim(2:NbLayers+1,:,:);
                            case 4
%                                 AllSubjects_Data(iROI).DATAPSC(:,:,SubjInd) = Data_ROI.MEDIANPSC(2:NbLayers+1.,:);
%                                 A=Data_ROI.LayerMedianPSC(2:NbLayers+1,:,:);
                        end
                    else
                        switch i
                            case 1
                                AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEAN(2:NbLayers+1,:);
                                A=Data_ROI.LayerMean(2:NbLayers+1,:,:);
                            case 2
                                AllSubjects_Data(iROI).DATARestrict(:,:,SubjInd) = Data_ROI.MEANRestrict(2:NbLayers+1,:);
                                A=Data_ROI.LayerMeanRestrict(2:NbLayers+1,:,:);
                            case 3
%                                 AllSubjects_Data(iROI).DATASkim(:,:,SubjInd) = Data_ROI.MEANSkim(2:NbLayers+1,:);
%                                 A=Data_ROI.LayerMeanSkim(2:NbLayers+1,:,:);
                            case 4
%                                 AllSubjects_Data(iROI).DATAPSC(:,:,SubjInd) = Data_ROI.MEANPSC(2:NbLayers+1.,:);
%                                 A=Data_ROI.LayerMeanPSC(2:NbLayers+1,:,:);
                        end
                    end
                    
                    % Compute main effect for each block
                    Blocks(:,:,1) = nanmean(A(:,:,[1 4]),3);
                    Blocks(:,:,2) = nanmean(A(:,:,[2 5]),3);
                    
                    % Compute main effect for each block
                    Blocks2 = A(:,:,1:3);
                    Blocks2(:,:,1) = A(:,:,3)-A(:,:,1) - (A(:,:,6)-A(:,:,4));
                    Blocks2(:,:,2) = A(:,:,3)-A(:,:,2) - (A(:,:,6)-A(:,:,5));
                    Blocks2(:,:,3) = A(:,:,2)-A(:,:,1) - (A(:,:,5)-A(:,:,4));
                    
                    switch i
                        case 1
                            AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = Blocks;
                            AllSubjects_Data(iROI).Interaction.Blocks.DATA{SubjInd} = Blocks2;
                        case 2
                            AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjInd} = Blocks;
                            AllSubjects_Data(iROI).InteractionRestrict.Blocks.DATA{SubjInd} = Blocks2;
                        case 3
%                             AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjInd} = Blocks;
%                             AllSubjects_Data(iROI).InteractionSkim.Blocks.DATA{SubjInd} = Blocks2;
                        case 4
%                             AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjInd} = Blocks;
%                             AllSubjects_Data(iROI).InteractionPSC.Blocks.DATA{SubjInd} = Blocks2;
                    end
                    
                    clear A Blocks Blocks2
                end
                
            end
            
            
        end
        
        clear Data_ROIs SubROI_Ind Data_ROI Data_ROI
        
    end
    
end

clear SubjInd ROI_Ind

cd(StartFolder)


%% Averages over subjects
for iROI=1:length(AllSubjects_Data)
    
    Include = logical(ones(NbSubject,1));
    % Include = ones(NbSubject,1); Include(4)=0; Include(11)=0; Include = logical(Include);
    % Include = all([~any(AllSubjects_Data(iROI).VoxelCount(:,2:end)<350,2) Include],2);
    % Include = [AllSubjects_Data(iROI).ROI_Coverage(:,1) ./ %AllSubjects_Data(iROI).ROI_Coverage(:,2)]>InclusionThreshold;
    
    AllSubjects_Data(iROI).Include=Include;
    
    for i=1:2
        switch i
            case 1
                A = AllSubjects_Data(iROI).DATA(:,:,Include);
            case 2
                A = AllSubjects_Data(iROI).DATARestrict(:,:,Include);
            case 3
%                 A = AllSubjects_Data(iROI).DATASkim(:,:,Include);
            case 4
%                 A = AllSubjects_Data(iROI).DATAPSC(:,:,Include);
        end
        % Compute main effects
        TEMP(:,1,:) = nanmean(A(:,[1 4],:),2);
        TEMP(:,2,:) = nanmean(A(:,[2 5],:),2);
        
        % Compute Interactions
        %   -AV-A(Att_V VS Att_A)
        %   -AV-V(Att_V VS Att_A)
        %   -A-V(Att_V VS Att_A)
        TEMP2(:,1,:) = A(:,3,:)-A(:,1,:) - (A(:,6,:)-A(:,4,:));
        TEMP2(:,2,:) = A(:,3,:)-A(:,2,:) - (A(:,6,:)-A(:,5,:));
        TEMP2(:,3,:) = A(:,2,:)-A(:,1,:) - (A(:,5,:)-A(:,4,:));
        
        switch i
            case 1
                AllSubjects_Data(iROI).MainEffects.DATA = TEMP;
                AllSubjects_Data(iROI).MainEffects.MEAN = nanmean(TEMP,3);
                AllSubjects_Data(iROI).MainEffects.STD = nanstd(TEMP,3);
                AllSubjects_Data(iROI).MainEffects.SEM = nansem(TEMP,3);
                
                AllSubjects_Data(iROI).Interaction.DATA = TEMP2;
                AllSubjects_Data(iROI).Interaction.MEAN = nanmean(TEMP2,3);
                AllSubjects_Data(iROI).Interaction.STD = nanstd(TEMP2,3);
                AllSubjects_Data(iROI).Interaction.SEM = nansem(TEMP2,3);
            case 2
                AllSubjects_Data(iROI).MainEffectsRestrict.DATA = TEMP;
                AllSubjects_Data(iROI).MainEffectsRestrict.MEAN = nanmean(TEMP,3);
                AllSubjects_Data(iROI).MainEffectsRestrict.STD = nanstd(TEMP,3);
                AllSubjects_Data(iROI).MainEffectsRestrict.SEM = nansem(TEMP,3);
                
                AllSubjects_Data(iROI).InteractionRestrict.DATA = TEMP2;
                AllSubjects_Data(iROI).InteractionRestrict.MEAN = nanmean(TEMP2,3);
                AllSubjects_Data(iROI).InteractionRestrict.STD = nanstd(TEMP2,3);
                AllSubjects_Data(iROI).InteractionRestrict.SEM = nansem(TEMP2,3);
            case 3
%                 AllSubjects_Data(iROI).MainEffectsSkim.DATA = TEMP;
%                 AllSubjects_Data(iROI).MainEffectsSkim.MEAN = nanmean(TEMP,3);
%                 AllSubjects_Data(iROI).MainEffectsSkim.STD = nanstd(TEMP,3);
%                 AllSubjects_Data(iROI).MainEffectsSkim.SEM = nansem(TEMP,3);
%                 
%                 AllSubjects_Data(iROI).InteractionSkim.DATA = TEMP2;
%                 AllSubjects_Data(iROI).InteractionSkim.MEAN = nanmean(TEMP2,3);
%                 AllSubjects_Data(iROI).InteractionSkim.STD = nanstd(TEMP2,3);
%                 AllSubjects_Data(iROI).InteractionSkim.SEM = nansem(TEMP2,3);
            case 4
%                 AllSubjects_Data(iROI).MainEffectsPSC.DATA = TEMP;
%                 AllSubjects_Data(iROI).MainEffectsPSC.MEAN = nanmean(TEMP,3);
%                 AllSubjects_Data(iROI).MainEffectsPSC.STD = nanstd(TEMP,3);
%                 AllSubjects_Data(iROI).MainEffectsPSC.SEM = nansem(TEMP,3);
                
%                 AllSubjects_Data(iROI).InteractionPSC.DATA = TEMP2;
%                 AllSubjects_Data(iROI).InteractionPSC.MEAN = nanmean(TEMP2,3);
%                 AllSubjects_Data(iROI).InteractionPSC.STD = nanstd(TEMP2,3);
%                 AllSubjects_Data(iROI).InteractionPSC.SEM = nansem(TEMP2,3);
        end
        
        
        clear A TEMP TEMP2
        
    end
    
end


%% Computes betas
for iROI=1:length(AllSubjects_Data)
    
    Name=AllSubjects_Data(iROI).name;
    
    fprintf('\nComputing betas for ROI %s\n', Name)
    
    SubjToInclude = find(AllSubjects_Data(iROI).Include);
    
    
    %% BOLD
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    AllSubjects_Data(iROI).Interaction.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    
    for i=1:2
        
        for SubjInd = 1:length(SubjToInclude)
            
            switch i
                case 1
                    Blocks = AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjToInclude(SubjInd)};
                case 2
                    Blocks = AllSubjects_Data(iROI).Interaction.Blocks.DATA{SubjToInclude(SubjInd)};
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
                            AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 2
                            AllSubjects_Data(iROI).Interaction.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                    end
                    
                    clear Y B
                end
                clear X Y B
                
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
    
    
    % Interaction average betas
    tmp=AllSubjects_Data(iROI).Interaction.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).Interaction.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).Interaction.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).Interaction.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).Interaction.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    clear tmp H P
    
    
    %% Restrict
    AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    
    for i=1:2
        
        for SubjInd = 1:length(SubjToInclude)
            
            switch i
                case 1
                    Blocks = AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjToInclude(SubjInd)};
                case 2
                    Blocks = AllSubjects_Data(iROI).InteractionRestrict.Blocks.DATA{SubjToInclude(SubjInd)};
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
                            AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 2
                            AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                    end
                    
                    clear Y B
                end
                clear X Y B
                
            end
            
        end
        
    end
    
    % Main Effects average betas
    tmp=AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Interaction average betas
    tmp=AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).InteractionRestrict.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    clear tmp H P
    
    
    %% Skim
%     AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     
%     for i=1:2
%         
%         for SubjInd = 1:length(SubjToInclude)
%             
%             switch i
%                 case 1
%                     Blocks = AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 2
%                     Blocks = AllSubjects_Data(iROI).InteractionSkim.Blocks.DATA{SubjToInclude(SubjInd)};
%             end
%             
%             if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
%                 
%                 for CondInd = 1:size(Blocks,3);
%                     
%                     Y = flipud(Blocks(:,:,CondInd));
%                     if any(isnan(Y(:)))
%                         [~,y]=find(isnan(Y));
%                         y=unique(y);
%                         Y(:,y)=[];
%                         clear y
%                     end
%                     
%                     if isempty(Y)
%                         B=nan(1,size(DesMat,2));
%                     else
%                         X=repmat(DesMat,size(Y,2),1);
%                         Y=Y(:);
%                         [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
%                     end
%                     
%                     switch i
%                         case 1
%                             AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 2
%                             AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                     end
%                     
%                     clear Y B
%                 end
%                 clear X Y B
%                 
%             end
%             
%         end
%         
%     end
%     
%     % Main Effects average betas
%     tmp=AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Interaction average betas
%     tmp=AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).InteractionSkim.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     clear tmp H P
    
    
    %% PSC
%     AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     
%     for i=1:2
%         
%         for SubjInd = 1:length(SubjToInclude)
%             
%             switch i
%                 case 1
%                     Blocks = AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 2
%                     Blocks = AllSubjects_Data(iROI).InteractionPSC.Blocks.DATA{SubjToInclude(SubjInd)};
%             end
%             
%             if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
%                 
%                 for CondInd = 1:size(Blocks,3);
%                     
%                     Y = flipud(Blocks(:,:,CondInd));
%                     if any(isnan(Y(:)))
%                         [~,y]=find(isnan(Y));
%                         y=unique(y);
%                         Y(:,y)=[];
%                         clear y
%                     end
%                     
%                     if isempty(Y)
%                         B=nan(1,size(DesMat,2));
%                     else
%                         X=repmat(DesMat,size(Y,2),1);
%                         Y=Y(:);
%                         [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
%                     end
%                     
%                     switch i
%                         case 1
%                             AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 2
%                             AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                     end
%                     
%                     clear Y B
%                 end
%                 clear X Y B
%                 
%             end
%             
%         end
%         
%     end
%     
%     % Main Effects average betas
%     tmp=AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Interaction average betas
%     tmp=AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).InteractionPSC.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     clear tmp H P
%     
end


%% Saves
fprintf('\nSaving\n')

cd(FigureFolder)
if Median
    if ANTs
        save(strcat('Data2_Median_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
    else
        save(strcat('Data2_Median_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
    end
else
    if ANTs
        save(strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
    else
        save(strcat('Data2_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
    end
end


cd(StartFolder)