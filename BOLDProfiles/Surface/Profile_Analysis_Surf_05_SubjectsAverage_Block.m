clear; close all; clc;

NbLayers = 6;
% NbLayers = NbLayers;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers+2) '_Layers'];

Median = 1;

InclusionThreshold = .5;


StartFolder = fullfile(pwd, '..', '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder, 'Figures', 'ProfilesSurface', strcat(num2str(NbLayers), '_layers'));
[~,~,~] = mkdir(FigureFolder);

% ROIs = {...
%     'A1';...
%     'A1_V_act';...
%     'A1_V_deact';...
%     'A1_A_act';...
%     'A1_A_deact';...
%     'PT';...
%     'PT_V_act';...
%     'PT_V_deact';...
%     'PT_A_act';...
%     'PT_A_deact';...
%     'V1';...
%     'V1_V_act';...
%     'V1_V_deact';...
%     'V1_A_act';...
%     'V1_A_deact';...
%     'V2-3';...
%     'V2-3_V_act';...
%     'V2-3_V_deact';...
%     'V2-3_A_act';...
%     'V2-3_A_deact';...
%     };

ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3';...
    'V1_act';...
    'V1_deact';...
    'V23_act';...
    'V23_deact';...
    };


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

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


%%
for iROI = 1:length(ROIs)
    AllSubjects_Data(iROI) = struct(...
        'name', ROIs{iROI}, ...
        'DATA', nan(NbLayers,6,NbSubject), ...
        'VertexCount', nan(NbSubject,NbLayers+1),...
        'Differential', struct('Blocks', struct('DATA', cell(1))));
end

%% Gets data for each subject
for SubjInd = 3 %1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Surfaces');
    
    
    for iROI=1%:size(ROIs,1)
        
        File2Load = fullfile(SubjectFolder,...
            strcat('Data_Surf_Block_', ROIs{iROI}, '_', num2str(NbLayers+2), '_layers.mat'));
        
        if ~exist(File2Load, 'file')
            error('File %s is missing.', File2Load)
        else
            
            load(File2Load, 'Data_ROI')
            
            if isfield(Data_ROI, 'MEAN')
                
                if Median
                    AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEDIAN(2:end-1,:);
                    AllSubjects_Data(iROI).NormDATA(:,:,SubjInd) = ...
                        Data_ROI.MEDIAN(2:end-1,:)./repmat(mean(Data_ROI.MEDIAN(2:end-1,:)), [size(Data_ROI.MEDIAN(2:end-1,:),1) 1 1]);
                else
                    AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEAN(2:end-1,:);
                    AllSubjects_Data(iROI).NormDATA(:,:,SubjInd) = ...
                        Data_ROI.MEAN(2:end-1,:)./repmat(mean(Data_ROI.MEAN(2:end-1,:)), [size(Data_ROI.MEAN(2:end-1,:),1) 1 1]);
                end
                
            end
            
            if isfield(Data_ROI, 'LayerMean')
                %% Compute main effect for each block
                if Median
                    A=Data_ROI.LayerMedian(2:end-1,:,:);
                else
                    A=Data_ROI.LayerMean(2:end-1,:,:); %#ok<*UNRCH>
                end
                
                Blocks = A(:,:,1:3);
                Blocks(:,:,end+1) = nanmean(A(:,:,1:3),3);
                Blocks(:,:,end+1:end+3) = A(:,:,4:6);
                Blocks(:,:,end+1) = nanmean(A(:,:,4:6),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[1 4]),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[2 5]),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[3 6]),3);
                Blocks(:,:,end+1) = nan(size(Blocks(:,:,end)));
                
                AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = Blocks;
                
                clear Blocks

                %% Compute differential for each block
                B = A;
                
                B(:,:,end+1:end+3) = B(:,:,4:6)-B(:,:,1:3);
                
                for j=1:3
                    TEMP(:,:,j) = B(:,:,j*3) - ( B(:,:,j*3-1)+B(:,:,j*3-2) ); %#ok<*SAGROW>
                end
                
                Blocks = B(:,:,1:3);
                Blocks(:,:,end+1) = TEMP(:,:,1);
                Blocks(:,:,end+1:end+3) = B(:,:,4:6);
                Blocks(:,:,end+1) = TEMP(:,:,2);
                Blocks(:,:,end+1:end+3) = B(:,:,7:9);
                Blocks(:,:,end+1) = TEMP(:,:,3);
                
                MainEffects{1} = mean(cat(3,Blocks(:,:,4),Blocks(:,:,8)),3);
                MainEffects{2} = mean(cat(3,Blocks(:,:,9),Blocks(:,:,10),Blocks(:,:,11)),3);
                
                AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjInd} = Blocks;
                AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                
                clear TEMP Blocks temp
                               
                
                %% Compute AV-V and AV-A for each block
                Blocks = A;
                Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,1), A(:,:,6)-A(:,:,4)),3);
                Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,2), A(:,:,6)-A(:,:,5)),3);
                
                AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjInd} = Blocks;
                
                clear TEMP Blocks temp
                
                
                %% Compute AV-V and AV-A for each block for each attention condition and also the interaction
                Blocks = A;
                Blocks(:,:,end+1) = A(:,:,3)-A(:,:,1);
                Blocks(:,:,end+1) = A(:,:,6)-A(:,:,4);
                Blocks(:,:,end+1) = (A(:,:,3)-A(:,:,1)) - (A(:,:,6)-A(:,:,4));
                Blocks(:,:,end+1) = A(:,:,3)-A(:,:,2);
                Blocks(:,:,end+1) = A(:,:,6)-A(:,:,5);
                Blocks(:,:,end+1) = (A(:,:,3)-A(:,:,2)) - (A(:,:,6)-A(:,:,5));
                
                AllSubjects_Data(iROI).BiVSUniSep.Blocks.DATA{SubjInd} = Blocks;
                
                clear TEMP Blocks temp
                
                
                
                %% Compute normalized main effect for each block
                A = A./repmat(mean(A), [size(A,1) 1 1]);
                
                Blocks = A(:,:,1:3);
                Blocks(:,:,end+1) = nanmean(A(:,:,1:3),3);
                Blocks(:,:,end+1:end+3) = A(:,:,4:6);
                Blocks(:,:,end+1) = nanmean(A(:,:,4:6),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[1 4]),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[2 5]),3);
                Blocks(:,:,end+1) = nanmean(A(:,:,[3 6]),3);
                Blocks(:,:,end+1) = nan(size(Blocks(:,:,end)));
                
                AllSubjects_Data(iROI).NormMainEffects.Blocks.DATA{SubjInd} = Blocks;
                
                clear Blocks
                
                                
                %% Compute normalized differential for each block
                B = A;
                
                B(:,:,end+1:end+3) = B(:,:,4:6)-B(:,:,1:3);
                
                for j=1:3
                    TEMP(:,:,j) = B(:,:,j*3) - ( B(:,:,j*3-1)+B(:,:,j*3-2) ); %#ok<*SAGROW>
                end
                
                Blocks = B(:,:,1:3);
                Blocks(:,:,end+1) = TEMP(:,:,1);
                Blocks(:,:,end+1:end+3) = B(:,:,4:6);
                Blocks(:,:,end+1) = TEMP(:,:,2);
                Blocks(:,:,end+1:end+3) = B(:,:,7:9);
                Blocks(:,:,end+1) = TEMP(:,:,3);
                
                MainEffects{1} = mean(cat(3,Blocks(:,:,4),Blocks(:,:,8)),3);
                MainEffects{2} = mean(cat(3,Blocks(:,:,9),Blocks(:,:,10),Blocks(:,:,11)),3);
                
                AllSubjects_Data(iROI).NormDifferential.Blocks.DATA{SubjInd} = Blocks;
                AllSubjects_Data(iROI).NormDifferential.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                
                clear TEMP Blocks temp
                
                
                %% Compute normalized AV-V and AV-A for each block
                Blocks = A;
                Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,1), A(:,:,6)-A(:,:,4)),3);
                Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,2), A(:,:,6)-A(:,:,5)),3);
                
                AllSubjects_Data(iROI).NormBiVSUni.Blocks.DATA{SubjInd} = Blocks;
                
                clear TEMP Blocks temp
                

            end
            
        end
        
        clear Data_ROIs SubROI_Ind Data_ROI
        
    end
    
end

clear SubjInd ROI_Ind

cd(StartFolder)


%% Averages over subjects
for iROI=1:length(AllSubjects_Data)
    
    Include = logical(ones(NbSubject,1)); %[AllSubjects_Data(iROI).ROI_Coverage(:,1) ./ ...
    %AllSubjects_Data(iROI).ROI_Coverage(:,2)]>InclusionThreshold;
    
    AllSubjects_Data(iROI).Include=Include;
    
    AllSubjects_Data(iROI).MEAN = nanmean(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).STD = nanstd(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).SEM = nansem(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    
    AllSubjects_Data(iROI).NormMEAN = nanmean(...
        AllSubjects_Data(iROI).NormDATA(:,:,Include),3);
    AllSubjects_Data(iROI).NormSTD = nanstd(...
        AllSubjects_Data(iROI).NormDATA(:,:,Include),3);
    AllSubjects_Data(iROI).NormSEM = nansem(...
        AllSubjects_Data(iROI).NormDATA(:,:,Include),3);
    
    %% BOLD
    A = AllSubjects_Data(iROI).DATA(:,:,Include);
    
    % Compute differential
    B=A;
    B = [B B(:,4:6,:)-B(:,1:3,:)]; %#ok<AGROW>
    
    for i=1:3
        TEMP(:,i,:) = B(:,i*3,:) - ( B(:,i*3-1,:)+B(:,i*3-2,:) );
    end
    B = [B(:,1:3,:) TEMP(:,1,:) B(:,4:6,:) TEMP(:,2,:) B(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    AllSubjects_Data(iROI).Differential.DATA = B;
    AllSubjects_Data(iROI).Differential.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).Differential.STD = nanstd(B,3);
    AllSubjects_Data(iROI).Differential.SEM = nansem(B,3);
    
    C(:,1,:) = nanmean(B(:,[4 8],:),2);
    C(:,2,:) = nanmean(B(:,9:11,:),2);
    
    AllSubjects_Data(iROI).Differential.MainEffect.DATA = C;
    AllSubjects_Data(iROI).Differential.MainEffect.MEAN = nanmean(C,3);
    AllSubjects_Data(iROI).Differential.MainEffect.STD = nanstd(C,3);
    AllSubjects_Data(iROI).Differential.MainEffect.SEM = nansem(C,3);
    
    clear B C
    
    % Compute main effects
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
    
    clear TEMP
    
    % Compute Bi - Uni
    TEMP = A(:,1:6,:);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
    
    AllSubjects_Data(iROI).BiVSUni.DATA = TEMP;
    AllSubjects_Data(iROI).BiVSUni.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).BiVSUni.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).BiVSUni.SEM = nansem(TEMP,3);
    
    clear TEMP
    
    
    % Compute Bi - Uni
    TEMP = A(:,1:6,:);
    TEMP(:,end+1,:) = A(:,3,:)-A(:,1,:);
    TEMP(:,end+1,:) = A(:,6,:)-A(:,4,:);
    TEMP(:,end+1,:) = (A(:,3,:)-A(:,1,:)) - (A(:,6,:)-A(:,4,:));
    TEMP(:,end+1,:) = A(:,3,:)-A(:,2,:);
    TEMP(:,end+1,:) = A(:,6,:)-A(:,5,:);
    TEMP(:,end+1,:) = (A(:,3,:)-A(:,2,:)) - (A(:,6,:)-A(:,5,:));
    
    AllSubjects_Data(iROI).BiVSUniSep.DATA = TEMP;
    AllSubjects_Data(iROI).BiVSUniSep.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).BiVSUniSep.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).BiVSUniSep.SEM = nansem(TEMP,3);
    
    clear A TEMP
    
    
    %% normalised BOLD
    A = AllSubjects_Data(iROI).NormDATA(:,:,Include);
    
    % Compute differential
    B=A;
    B = [B B(:,4:6,:)-B(:,1:3,:)]; %#ok<AGROW>
    
    for i=1:3
        TEMP(:,i,:) = B(:,i*3,:) - ( B(:,i*3-1,:)+B(:,i*3-2,:) );
    end
    B = [B(:,1:3,:) TEMP(:,1,:) B(:,4:6,:) TEMP(:,2,:) B(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    AllSubjects_Data(iROI).NormDifferential.DATA = B;
    AllSubjects_Data(iROI).NormDifferential.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).NormDifferential.STD = nanstd(B,3);
    AllSubjects_Data(iROI).NormDifferential.SEM = nansem(B,3);
    
    C(:,1,:) = nanmean(B(:,[4 8],:),2);
    C(:,2,:) = nanmean(B(:,9:11,:),2);
    
    AllSubjects_Data(iROI).NormDifferential.MainEffect.DATA = C;
    AllSubjects_Data(iROI).NormDifferential.MainEffect.MEAN = nanmean(C,3);
    AllSubjects_Data(iROI).NormDifferential.MainEffect.STD = nanstd(C,3);
    AllSubjects_Data(iROI).NormDifferential.MainEffect.SEM = nansem(C,3);
    
    clear B C
    
    % Compute main effects
    TEMP(:,1:3,:) = A(:,1:3,:);
    TEMP(:,end+1,:) = nanmean(A(:,1:3,:),2);
    TEMP(:,end+1:end+3,:) = A(:,4:6,:);
    TEMP(:,end+1,:) = nanmean(A(:,4:6,:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[1 4],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[2 5],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[3 6],:),2);
    TEMP(:,end+1,:) = nan(size(TEMP(:,end,:)));
    
    AllSubjects_Data(iROI).NormMainEffects.DATA = TEMP;
    AllSubjects_Data(iROI).NormMainEffects.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).NormMainEffects.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).NormMainEffects.SEM = nansem(TEMP,3);
    
    clear TEMP
    
    % Compute Bi - Uni
    TEMP = A(:,1:6,:);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
    
    AllSubjects_Data(iROI).NormBiVSUni.DATA = TEMP;
    AllSubjects_Data(iROI).NormBiVSUni.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).NormBiVSUni.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).NormBiVSUni.SEM = nansem(TEMP,3);
    
    clear A TEMP
    
end

%% Computes betas
for iROI=1:length(AllSubjects_Data)
    
    
    Name=AllSubjects_Data(iROI).name;
    
    fprintf('\nComputing betas for ROI %s\n', Name)
    
    SubjToInclude = find(AllSubjects_Data(iROI).Include);
    
    
    %% BOLD
    AllSubjects_Data(iROI).Differential.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA = ...
        nan(size(DesMat,2),2,size(SubjectList,1));
    AllSubjects_Data(iROI).MainEffects.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    
    for i=1:4
        
        for SubjInd = 1:length(SubjToInclude)
            
            switch i
                case 1
                    Blocks = AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjToInclude(SubjInd)};
                case 2
                    Blocks = AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjToInclude(SubjInd)};
                case 3
                    Blocks = AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjToInclude(SubjInd)};
                case 4
                    Blocks = AllSubjects_Data(iROI).BiVSUniSep.Blocks.DATA{SubjToInclude(SubjInd)};                    
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
                        case 4
                            AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;                            
                    end
                    
                    clear Y B
                end
                clear X Y B
                
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
    
    
    % Bi VS Uni average betas
    tmp=AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).BiVSUniSep.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end    
    
    
    clear tmp H P
end


%% Saves
fprintf('\nSaving\n')

cd(FigureFolder)
if Median
    save(strcat('Data_Surf_Median_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
else
    save(strcat('Data_Surf_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
end
cd(StartFolder)