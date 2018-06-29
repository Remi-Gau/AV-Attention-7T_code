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
                
                AllSubjects_Data(iROI).DATA(:,:,SubjInd) = nan(NbLayers,6);
                AllSubjects_Data(iROI).DATARestrict(:,:,SubjInd) = nan(NbLayers,6);
                AllSubjects_Data(iROI).DATASkim(:,:,SubjInd) = nan(NbLayers,6);
                AllSubjects_Data(iROI).DATAPSC(:,:,SubjInd) = nan(NbLayers,6);
                
                AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
%                 AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                
                AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjInd} = {nan(NbLayers,12) nan(NbLayers,12)};
                AllSubjects_Data(iROI).DifferentialRestrict.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.DATA{SubjInd} = {nan(NbLayers,12) nan(NbLayers,12)};
                AllSubjects_Data(iROI).DifferentialSkim.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
                AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.DATA{SubjInd} = {nan(NbLayers,12) nan(NbLayers,12)};
%                 AllSubjects_Data(iROI).DifferentialPSC.Blocks.DATA{SubjInd} = nan(NbLayers,12,12);
%                 AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.DATA{SubjInd} = {nan(NbLayers,12) nan(NbLayers,12)};
                
                AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjInd} = nan(NbLayers,12,8);
                AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.DATA{SubjInd} = nan(NbLayers,12,8);
                AllSubjects_Data(iROI).BiVSUniSkim.Blocks.DATA{SubjInd} = nan(NbLayers,12,8);
%                 AllSubjects_Data(iROI).BiVSUniPSC.Blocks.DATA{SubjInd} = nan(NbLayers,12,8);
                
            else
                
                AllSubjects_Data(iROI).VoxelCount(SubjInd,1:NbLayers+1) = Data_ROI.Voxelcount;
                
                
                if Median
                    AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEDIAN(2:NbLayers+1,:);
                    AllSubjects_Data(iROI).DATARestrict(:,:,SubjInd) = Data_ROI.MEDIANRestrict(2:NbLayers+1,:);
%                     AllSubjects_Data(iROI).DATASkim(:,:,SubjInd) = Data_ROI.MEDIANSkim(2:NbLayers+1,:);
%                     AllSubjects_Data(iROI).DATAPSC(:,:,SubjInd) = Data_ROI.MEDIANPSC(2:NbLayers+1,:);
                else
                    AllSubjects_Data(iROI).DATA(:,:,SubjInd) = Data_ROI.MEAN(2:NbLayers+1,:);
                    AllSubjects_Data(iROI).DATARestrict(:,:,SubjInd) = Data_ROI.MEANRestrict(2:NbLayers+1,:);
%                     AllSubjects_Data(iROI).DATASkim(:,:,SubjInd) = Data_ROI.MEANSkim(2:NbLayers+1,:);
%                     AllSubjects_Data(iROI).DATAPSC(:,:,SubjInd) = Data_ROI.MEANPSC(2:NbLayers+1,:);
                end
                
                for i=1:3
                    
                    %% Compute main effect for each block
                    if Median
                        if i==1
                            A=Data_ROI.LayerMedian(2:NbLayers+1,:,:);
                        elseif i==2
                            A=Data_ROI.LayerMedianRestrict(2:NbLayers+1,:,:);
                        elseif i==3
%                             A=Data_ROI.LayerMedianSkim(2:NbLayers+1,:,:);
                        else
%                             A=Data_ROI.LayerMedianPSC(2:NbLayers+1,:,:);
                        end
                    else
                        if i==1
                            A=Data_ROI.LayerMean(2:NbLayers+1,:,:);
                        elseif i==2
                            A=Data_ROI.LayerMeanRestrict(2:NbLayers+1,:,:);
                        elseif i==3
%                             A=Data_ROI.LayerMeanSkim(2:NbLayers+1,:,:);
                        else
%                             A=Data_ROI.LayerMeanPSC(2:NbLayers+1,:,:);
                        end
                    end
                    
                    Blocks = A(:,:,1:3);
                    Blocks(:,:,end+1) = nanmean(A(:,:,1:3),3);
                    Blocks(:,:,end+1:end+3) = A(:,:,4:6);
                    Blocks(:,:,end+1) = nanmean(A(:,:,4:6),3);
                    Blocks(:,:,end+1) = nanmean(A(:,:,[1 4]),3);
                    Blocks(:,:,end+1) = nanmean(A(:,:,[2 5]),3);
                    Blocks(:,:,end+1) = nanmean(A(:,:,[3 6]),3);
                    Blocks(:,:,end+1) = nan(size(Blocks(:,:,end)));
                    
                    if i==1
                        AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = Blocks;
                    elseif i==2
                        AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjInd} = Blocks;
                    elseif i==3
%                         AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjInd} = Blocks;
                    else
%                         AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjInd} = Blocks;
                    end
                    
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
                    
                    if i==1
                        AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjInd} = Blocks;
                        AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                    elseif i==2
                        AllSubjects_Data(iROI).DifferentialRestrict.Blocks.DATA{SubjInd} = Blocks;
                        AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                    elseif i==3
                        AllSubjects_Data(iROI).DifferentialSkim.Blocks.DATA{SubjInd} = Blocks;
                        AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                    else
%                         AllSubjects_Data(iROI).DifferentialPSC.Blocks.DATA{SubjInd} = Blocks;
%                         AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
                    end
                    
                    clear TEMP Blocks temp
                    
                    
                    %% Compute AV-V and AV-A for each block
                    Blocks = A;
                    Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,1), A(:,:,6)-A(:,:,4)),3);
                    Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,2), A(:,:,6)-A(:,:,5)),3);
                    
                    if i==1
                        AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjInd} = Blocks;
                    elseif i==2
                        AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.DATA{SubjInd} = Blocks;
                    elseif i==3
                        AllSubjects_Data(iROI).BiVSUniSkim.Blocks.DATA{SubjInd} = Blocks;
                    else
%                         AllSubjects_Data(iROI).BiVSUniPSC.Blocks.DATA{SubjInd} = Blocks;
                    end
                    
                    clear TEMP Blocks temp
                    
                    
                end
                
            end
        end
        
        clear Data_ROIs SubROI_Ind Data_ROI
        
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
    
    AllSubjects_Data(iROI).MEAN = nanmean(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).STD = nanstd(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    AllSubjects_Data(iROI).SEM = nansem(...
        AllSubjects_Data(iROI).DATA(:,:,Include),3);
    
    AllSubjects_Data(iROI).MEANRestrict = nanmean(...
        AllSubjects_Data(iROI).DATARestrict(:,:,Include),3);
    AllSubjects_Data(iROI).STDRestrict = nanstd(...
        AllSubjects_Data(iROI).DATARestrict(:,:,Include),3);
    AllSubjects_Data(iROI).SEMRestrict = nansem(...
        AllSubjects_Data(iROI).DATARestrict(:,:,Include),3);
    
%     AllSubjects_Data(iROI).MEANSkim = nanmean(...
%         AllSubjects_Data(iROI).DATASkim(:,:,Include),3);
%     AllSubjects_Data(iROI).STDSkim = nanstd(...
%         AllSubjects_Data(iROI).DATASkim(:,:,Include),3);
%     AllSubjects_Data(iROI).SEMSkim = nansem(...
%         AllSubjects_Data(iROI).DATASkim(:,:,Include),3);
    
%     AllSubjects_Data(iROI).MEANPSC = nanmean(...
%         AllSubjects_Data(iROI).DATAPSC(:,:,Include),3);
%     AllSubjects_Data(iROI).STDPSC = nanstd(...
%         AllSubjects_Data(iROI).DATAPSC(:,:,Include),3);
%     AllSubjects_Data(iROI).SEMPSC = nansem(...
%         AllSubjects_Data(iROI).DATAPSC(:,:,Include),3);
    
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
    
    % Compute derivative
    B = diff(A,1,1);
    AllSubjects_Data(iROI).Derivative.DATA = B;
    AllSubjects_Data(iROI).Derivative.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).Derivative.STD = nanstd(B,3);
    AllSubjects_Data(iROI).Derivative.SEM = nansem(B,3);
    clear B
    
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
    
    clear A TEMP
    
    
    %% Restrict
    A = AllSubjects_Data(iROI).DATARestrict(:,:,Include);
    
    % Compute differential
    B=A;
    B = [B B(:,4:6,:)-B(:,1:3,:)]; %#ok<AGROW>
    
    for i=1:3
        TEMP(:,i,:) = B(:,i*3,:) - ( B(:,i*3-1,:)+B(:,i*3-2,:) );
    end
    B = [B(:,1:3,:) TEMP(:,1,:) B(:,4:6,:) TEMP(:,2,:) B(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    AllSubjects_Data(iROI).DifferentialRestrict.DATA = B;
    AllSubjects_Data(iROI).DifferentialRestrict.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).DifferentialRestrict.STD = nanstd(B,3);
    AllSubjects_Data(iROI).DifferentialRestrict.SEM = nansem(B,3);
    
    C(:,1,:) = nanmean(B(:,[4 8],:),2);
    C(:,2,:) = nanmean(B(:,9:11,:),2);
    
    AllSubjects_Data(iROI).DifferentialRestrict.MainEffect.DATA = C;
    AllSubjects_Data(iROI).DifferentialRestrict.MainEffect.MEAN = nanmean(C,3);
    AllSubjects_Data(iROI).DifferentialRestrict.MainEffect.STD = nanstd(C,3);
    AllSubjects_Data(iROI).DifferentialRestrict.MainEffect.SEM = nansem(C,3);
    
    clear B C
    
    % Compute derivative
    B = diff(A,1,1);
    AllSubjects_Data(iROI).DerivativeRestrict.DATA = B;
    AllSubjects_Data(iROI).DerivativeRestrict.MEAN = nanmean(B,3);
    AllSubjects_Data(iROI).DerivativeRestrict.STD = nanstd(B,3);
    AllSubjects_Data(iROI).DerivativeRestrict.SEM = nansem(B,3);
    clear B
    
    % Compute main effects
    TEMP(:,1:3,:) = A(:,1:3,:);
    TEMP(:,end+1,:) = nanmean(A(:,1:3,:),2);
    TEMP(:,end+1:end+3,:) = A(:,4:6,:);
    TEMP(:,end+1,:) = nanmean(A(:,4:6,:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[1 4],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[2 5],:),2);
    TEMP(:,end+1,:) = nanmean(A(:,[3 6],:),2);
    TEMP(:,end+1,:) = nan(size(TEMP(:,end,:)));
    
    AllSubjects_Data(iROI).MainEffectsRestrict.DATA = TEMP;
    AllSubjects_Data(iROI).MainEffectsRestrict.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).MainEffectsRestrict.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).MainEffectsRestrict.SEM = nansem(TEMP,3);
    
    clear TEMP
    
    % Compute Bi - Uni
    TEMP = A(:,1:6,:);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
    TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
    
    AllSubjects_Data(iROI).BiVSUniRestrict.DATA = TEMP;
    AllSubjects_Data(iROI).BiVSUniRestrict.MEAN = nanmean(TEMP,3);
    AllSubjects_Data(iROI).BiVSUniRestrict.STD = nanstd(TEMP,3);
    AllSubjects_Data(iROI).BiVSUniRestrict.SEM = nansem(TEMP,3);
    
    clear A TEMP
    
    
    %% Skim
%     A = AllSubjects_Data(iROI).DATASkim(:,:,Include);
%     
%     % Compute differential
%     B=A;
%     B = [B B(:,4:6,:)-B(:,1:3,:)]; %#ok<AGROW>
%     
%     for i=1:3
%         TEMP(:,i,:) = B(:,i*3,:) - ( B(:,i*3-1,:)+B(:,i*3-2,:) );
%     end
%     B = [B(:,1:3,:) TEMP(:,1,:) B(:,4:6,:) TEMP(:,2,:) B(:,7:9,:) TEMP(:,3,:)];
%     clear TEMP
%     
%     AllSubjects_Data(iROI).DifferentialSkim.DATA = B;
%     AllSubjects_Data(iROI).DifferentialSkim.MEAN = nanmean(B,3);
%     AllSubjects_Data(iROI).DifferentialSkim.STD = nanstd(B,3);
%     AllSubjects_Data(iROI).DifferentialSkim.SEM = nansem(B,3);
%     
%     C(:,1,:) = nanmean(B(:,[4 8],:),2);
%     C(:,2,:) = nanmean(B(:,9:11,:),2);
%     
%     AllSubjects_Data(iROI).DifferentialSkim.MainEffect.DATA = C;
%     AllSubjects_Data(iROI).DifferentialSkim.MainEffect.MEAN = nanmean(C,3);
%     AllSubjects_Data(iROI).DifferentialSkim.MainEffect.STD = nanstd(C,3);
%     AllSubjects_Data(iROI).DifferentialSkim.MainEffect.SEM = nansem(C,3);
%     
%     clear B C
%     
%     % Compute derivative
%     B = diff(A,1,1);
%     AllSubjects_Data(iROI).DerivativeSkim.DATA = B;
%     AllSubjects_Data(iROI).DerivativeSkim.MEAN = nanmean(B,3);
%     AllSubjects_Data(iROI).DerivativeSkim.STD = nanstd(B,3);
%     AllSubjects_Data(iROI).DerivativeSkim.SEM = nansem(B,3);
%     clear B
%     
%     % Compute main effects
%     TEMP(:,1:3,:) = A(:,1:3,:);
%     TEMP(:,end+1,:) = nanmean(A(:,1:3,:),2);
%     TEMP(:,end+1:end+3,:) = A(:,4:6,:);
%     TEMP(:,end+1,:) = nanmean(A(:,4:6,:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[1 4],:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[2 5],:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[3 6],:),2);
%     TEMP(:,end+1,:) = nan(size(TEMP(:,end,:)));
%     
%     AllSubjects_Data(iROI).MainEffectsSkim.DATA = TEMP;
%     AllSubjects_Data(iROI).MainEffectsSkim.MEAN = nanmean(TEMP,3);
%     AllSubjects_Data(iROI).MainEffectsSkim.STD = nanstd(TEMP,3);
%     AllSubjects_Data(iROI).MainEffectsSkim.SEM = nansem(TEMP,3);
%     
%     clear TEMP
%     
%     % Compute Bi - Uni
%     TEMP = A(:,1:6,:);
%     TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
%     TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
%     
%     AllSubjects_Data(iROI).BiVSUniSkim.DATA = TEMP;
%     AllSubjects_Data(iROI).BiVSUniSkim.MEAN = nanmean(TEMP,3);
%     AllSubjects_Data(iROI).BiVSUniSkim.STD = nanstd(TEMP,3);
%     AllSubjects_Data(iROI).BiVSUniSkim.SEM = nansem(TEMP,3);
%     
%     clear A TEMP
    
    
    %% PSC
%     A = AllSubjects_Data(iROI).DATAPSC(:,:,Include);
%     
%     % Compute differential
%     B=A;
%     B = [B B(:,4:6,:)-B(:,1:3,:)]; %#ok<AGROW>
%     
%     for i=1:3
%         TEMP(:,i,:) = B(:,i*3,:) - ( B(:,i*3-1,:)+B(:,i*3-2,:) );
%     end
%     B = [B(:,1:3,:) TEMP(:,1,:) B(:,4:6,:) TEMP(:,2,:) B(:,7:9,:) TEMP(:,3,:)];
%     clear TEMP
%     
%     AllSubjects_Data(iROI).DifferentialPSC.DATA = B;
%     AllSubjects_Data(iROI).DifferentialPSC.MEAN = nanmean(B,3);
%     AllSubjects_Data(iROI).DifferentialPSC.STD = nanstd(B,3);
%     AllSubjects_Data(iROI).DifferentialPSC.SEM = nansem(B,3);
%     
%     C(:,1,:) = nanmean(B(:,[4 8],:),2);
%     C(:,2,:) = nanmean(B(:,9:11,:),2);
%     
%     AllSubjects_Data(iROI).DifferentialPSC.MainEffect.DATA = C;
%     AllSubjects_Data(iROI).DifferentialPSC.MainEffect.MEAN = nanmean(C,3);
%     AllSubjects_Data(iROI).DifferentialPSC.MainEffect.STD = nanstd(C,3);
%     AllSubjects_Data(iROI).DifferentialPSC.MainEffect.SEM = nansem(C,3);
%     
%     clear B C
%     
%     % Compute derivative
%     B = diff(A,1,1);
%     AllSubjects_Data(iROI).DerivativePSC.DATA = B;
%     AllSubjects_Data(iROI).DerivativePSC.MEAN = nanmean(B,3);
%     AllSubjects_Data(iROI).DerivativePSC.STD = nanstd(B,3);
%     AllSubjects_Data(iROI).DerivativePSC.SEM = nansem(B,3);
%     clear B
%     
%     % Compute main effects
%     TEMP(:,1:3,:) = A(:,1:3,:);
%     TEMP(:,end+1,:) = nanmean(A(:,1:3,:),2);
%     TEMP(:,end+1:end+3,:) = A(:,4:6,:);
%     TEMP(:,end+1,:) = nanmean(A(:,4:6,:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[1 4],:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[2 5],:),2);
%     TEMP(:,end+1,:) = nanmean(A(:,[3 6],:),2);
%     TEMP(:,end+1,:) = nan(size(TEMP(:,end,:)));
%     
%     AllSubjects_Data(iROI).MainEffectsPSC.DATA = TEMP;
%     AllSubjects_Data(iROI).MainEffectsPSC.MEAN = nanmean(TEMP,3);
%     AllSubjects_Data(iROI).MainEffectsPSC.STD = nanstd(TEMP,3);
%     AllSubjects_Data(iROI).MainEffectsPSC.SEM = nansem(TEMP,3);
%     
%     clear TEMP
%     
%     % Compute Bi - Uni
%     TEMP = A(:,1:6,:);
%     TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),4);
%     TEMP(:,end+1,:) = mean(cat(4,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),4);
%     
%     AllSubjects_Data(iROI).BiVSUniPSC.DATA = TEMP;
%     AllSubjects_Data(iROI).BiVSUniPSC.MEAN = nanmean(TEMP,3);
%     AllSubjects_Data(iROI).BiVSUniPSC.STD = nanstd(TEMP,3);
%     AllSubjects_Data(iROI).BiVSUniPSC.SEM = nansem(TEMP,3);
%     
%     clear A TEMP
    
    
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
    
    
    clear tmp H P
    
    
    %% Restrict
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.DATA = ...
        nan(size(DesMat,2),2,size(SubjectList,1));
    AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA = ...
        nan(size(DesMat,2),12,size(SubjectList,1));
    
    
    for i=1:3
        
        for SubjInd = 1:length(SubjToInclude)
            
            switch i
                case 1
                    Blocks = AllSubjects_Data(iROI).DifferentialRestrict.Blocks.DATA{SubjToInclude(SubjInd)};
                case 2
                    Blocks = AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.DATA{SubjToInclude(SubjInd)};
                case 3
                    Blocks = AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.DATA{SubjToInclude(SubjInd)};
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
                            AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 2
                            AllSubjects_Data(iROI).MainEffectsRestrict.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                        case 3
                            AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
                    end
                    
                    clear Y B
                end
                clear X Y B
                
            end
            
            % Betas for main effect of differential
            if i==1
                for j=1:numel(AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)})
                    Blocks = AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)}{j};
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
                        
                        AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.DATA(:,j,SubjToInclude(SubjInd))=B;
                        
                        clear Y B
                    end
                    clear X
                end
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
    
    
    % Differentials average betas
    tmp=AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).DifferentialRestrict.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Differentials main effects average betas
    tmp=AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.DATA;
    
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).DifferentialRestrict.Blocks.MainEffects.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    % Bi VS Uni average betas
    tmp=AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.DATA;
    
    AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.MEAN=nanmean(tmp, 3);
    AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.STD=nanstd(tmp, 3);
    AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.SEM=nansem(tmp, 3);
    
    % T-Test
    for CondInd = 1:size(tmp,2)
        for BetaInd=1:size(tmp,1)
            [H,P] = ttest(tmp(BetaInd,CondInd,:));
            AllSubjects_Data(iROI).BiVSUniRestrict.Blocks.Beta.P(BetaInd,CondInd)=P;
        end
    end
    
    
    clear tmp H P
    
    
    %% Skim
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.DATA = ...
%         nan(size(DesMat,2),2,size(SubjectList,1));
%     AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     
%     for i=1:3
%         
%         for SubjInd = 1:length(SubjToInclude)
%             
%             switch i
%                 case 1
%                     Blocks = AllSubjects_Data(iROI).DifferentialSkim.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 2
%                     Blocks = AllSubjects_Data(iROI).MainEffectsSkim.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 3
%                     Blocks = AllSubjects_Data(iROI).BiVSUniSkim.Blocks.DATA{SubjToInclude(SubjInd)};
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
%                             AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 2
%                             AllSubjects_Data(iROI).MainEffectsSkim.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 3
%                             AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                     end
%                     
%                     clear Y B
%                 end
%                 clear X Y B
%                 
%             end
%             
%             % Betas for main effect of differential
%             if i==1
%                 for j=1:numel(AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)})
%                     Blocks = AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)}{j};
%                     if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
%                         Y = flipud(Blocks(:,:));
%                         if any(isnan(Y(:)))
%                             [~,y]=find(isnan(Y));
%                             y=unique(y);
%                             Y(:,y)=[];
%                             clear y
%                         end
%                         
%                         if isempty(Y)
%                             B=nan(1,size(DesMat,2));
%                         else
%                             X=repmat(DesMat,size(Y,2),1);
%                             Y=Y(:);
%                             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
%                         end
%                         
%                         AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.DATA(:,j,SubjToInclude(SubjInd))=B;
%                         
%                         clear Y B
%                     end
%                     clear X
%                 end
%             end
%             
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
%     % Differentials average betas
%     tmp=AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).DifferentialSkim.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Differentials main effects average betas
%     tmp=AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.DATA;
%     
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).DifferentialSkim.Blocks.MainEffects.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Bi VS Uni average betas
%     tmp=AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).BiVSUniSkim.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     clear tmp H P
    
    
    
    %% PSC
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.DATA = ...
%         nan(size(DesMat,2),2,size(SubjectList,1));
%     AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.DATA = ...
%         nan(size(DesMat,2),12,size(SubjectList,1));
%     
%     
%     for i=1:3
%         
%         for SubjInd = 1:length(SubjToInclude)
%             
%             switch i
%                 case 1
%                     Blocks = AllSubjects_Data(iROI).DifferentialPSC.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 2
%                     Blocks = AllSubjects_Data(iROI).MainEffectsPSC.Blocks.DATA{SubjToInclude(SubjInd)};
%                 case 3
%                     Blocks = AllSubjects_Data(iROI).BiVSUniPSC.Blocks.DATA{SubjToInclude(SubjInd)};
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
%                             AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 2
%                             AllSubjects_Data(iROI).MainEffectsPSC.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                         case 3
%                             AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.DATA(:,CondInd,SubjToInclude(SubjInd))=B;
%                     end
%                     
%                     clear Y B
%                 end
%                 clear X Y B
%                 
%             end
%             
%             % Betas for main effect of differential
%             if i==1
%                 for j=1:numel(AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)})
%                     Blocks = AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.DATA{SubjToInclude(SubjInd)}{j};
%                     if ~all(isnan(Blocks(:))) || ~isempty(Blocks)
%                         Y = flipud(Blocks(:,:));
%                         if any(isnan(Y(:)))
%                             [~,y]=find(isnan(Y));
%                             y=unique(y);
%                             Y(:,y)=[];
%                             clear y
%                         end
%                         
%                         if isempty(Y)
%                             B=nan(1,size(DesMat,2));
%                         else
%                             X=repmat(DesMat,size(Y,2),1);
%                             Y=Y(:);
%                             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
%                         end
%                         
%                         AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.DATA(:,j,SubjToInclude(SubjInd))=B;
%                         
%                         clear Y B
%                     end
%                     clear X
%                 end
%             end
%             
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
%     % Differentials average betas
%     tmp=AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).DifferentialPSC.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Differentials main effects average betas
%     tmp=AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.DATA;
%     
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).DifferentialPSC.Blocks.MainEffects.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     % Bi VS Uni average betas
%     tmp=AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.DATA;
%     
%     AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.MEAN=nanmean(tmp, 3);
%     AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.STD=nanstd(tmp, 3);
%     AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.SEM=nansem(tmp, 3);
%     
%     % T-Test
%     for CondInd = 1:size(tmp,2)
%         for BetaInd=1:size(tmp,1)
%             [H,P] = ttest(tmp(BetaInd,CondInd,:));
%             AllSubjects_Data(iROI).BiVSUniPSC.Blocks.Beta.P(BetaInd,CondInd)=P;
%         end
%     end
%     
%     
%     clear tmp H P
    
    
end




%% Saves
fprintf('\nSaving\n')

cd(FigureFolder)
if Median
    if ANTs
        save(strcat('Data_Median_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
    else
        save(strcat('Data_Median_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
    end
else
    if ANTs
        save(strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '_ANTs.mat'))
    else
        save(strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
        % save(strcat('Data_Block_QuadGLM_', num2str(NbLayers), '_Layers', '.mat'))
    end
end
cd(StartFolder)
