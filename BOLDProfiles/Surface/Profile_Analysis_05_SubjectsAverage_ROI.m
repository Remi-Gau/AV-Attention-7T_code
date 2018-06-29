clear; close all; clc;

NbLayers = 6;
% NbLayers = NbLayers;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers+2) '_Layers'];

Median = 1;

InclusionThreshold = .5;


StartFolder = fullfile(pwd, '..', '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

%%
for iROI = 1:length(ROIs)
    AllSubjects_Data(iROI) = struct(...
        'name', ROIs{iROI}, ...
        'DATA', nan(NbSubject,6), ...
        'Differential', struct('Blocks', struct('DATA', cell(1))));
end

%% Gets data for each subject
for SubjInd = 1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Surfaces');
    
    
    for iROI=1:size(ROIs,1)
        
        File2Load = fullfile(SubjectFolder,...
            strcat('Data_Surf_Block_', ROIs{iROI}, '_', num2str(NbLayers+2), '_layers.mat'));
        
        if ~exist(File2Load, 'file')
            error('File %s is missing.', File2Load)
        else
            
            load(File2Load, 'Data_ROI')
            
            
            AllSubjects_Data(iROI).DATA(SubjInd,:) = mean(Data_ROI.MEDIAN(2:7,:));
            
            
            %% Compute main effect for each block
            A=mean(Data_ROI.LayerMedian(2:7,:,:),1);
            
            
            Blocks = A(:,:,1:3);
            Blocks(:,:,end+1) = nanmean(A(:,:,1:3),3);
            Blocks(:,:,end+1:end+3) = A(:,:,4:6);
            Blocks(:,:,end+1) = nanmean(A(:,:,4:6),3);
            Blocks(:,:,end+1) = nanmean(A(:,:,[1 4]),3);
            Blocks(:,:,end+1) = nanmean(A(:,:,[2 5]),3);
            Blocks(:,:,end+1) = nanmean(A(:,:,[3 6]),3);
            Blocks(:,:,end+1) = nan(size(Blocks(:,:,end)));
            
            AllSubjects_Data(iROI).MainEffects.Blocks.DATA{SubjInd} = squeeze(Blocks)';
            
            
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
            
            
            AllSubjects_Data(iROI).Differential.Blocks.DATA{SubjInd} = squeeze(Blocks)';
            AllSubjects_Data(iROI).Differential.Blocks.MainEffects.DATA{SubjInd} = MainEffects;
            
            
            clear TEMP Blocks temp
            
            
            %% Compute AV-V and AV-A for each block
            Blocks = A;
            Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,1), A(:,:,6)-A(:,:,4)),3);
            Blocks(:,:,end+1) = mean(cat(3,A(:,:,3)-A(:,:,2), A(:,:,6)-A(:,:,5)),3);
            
            
            AllSubjects_Data(iROI).BiVSUni.Blocks.DATA{SubjInd} = squeeze(Blocks)';
            
            
            %% Compute Interaction for each block
            Blocks2 = A(:,:,1:3);
            Blocks2(:,:,1) = A(:,:,3)-A(:,:,1) - (A(:,:,6)-A(:,:,4));
            Blocks2(:,:,2) = A(:,:,3)-A(:,:,2) - (A(:,:,6)-A(:,:,5));
            Blocks2(:,:,3) = A(:,:,2)-A(:,:,1) - (A(:,:,5)-A(:,:,4));
            
            AllSubjects_Data(iROI).Interaction.Blocks.DATA{SubjInd} = squeeze(Blocks2)';
            
            
            clear TEMP Blocks temp

        end
        
        clear Data_ROIs SubROI_Ind Data_ROI
        
    end
    
end

clear SubjInd ROI_Ind

cd(StartFolder)


%% Averages over subjects
for iROI=1:length(AllSubjects_Data)
    
    Include = logical(ones(NbSubject,1));
    
    AllSubjects_Data(iROI).Include=Include;
    
    AllSubjects_Data(iROI).MEAN = nanmean(AllSubjects_Data(iROI).DATA(Include,:));
    AllSubjects_Data(iROI).STD = nanstd(AllSubjects_Data(iROI).DATA(Include,:));
    AllSubjects_Data(iROI).SEM = nansem(AllSubjects_Data(iROI).DATA(Include,:));
    
    
    %% BOLD
    A = AllSubjects_Data(iROI).DATA(Include,:);
    
    % Compute differential
    B=A;
    B = [B B(:,4:6)-B(:,1:3)]; %#ok<AGROW>
    
    for i=1:3
        TEMP(:,i) = B(:,i*3) - ( B(:,i*3-1)+B(:,i*3-2) );
    end
    B = [B(:,1:3) TEMP(:,1) B(:,4:6) TEMP(:,2) B(:,7:9) TEMP(:,3)];
    clear TEMP
    
    AllSubjects_Data(iROI).Differential.DATA = B;
    AllSubjects_Data(iROI).Differential.MEAN = nanmean(B);
    AllSubjects_Data(iROI).Differential.STD = nanstd(B);
    AllSubjects_Data(iROI).Differential.SEM = nansem(B);
    
    C(:,1) = nanmean(B(:,[4 8]),2);
    C(:,2) = nanmean(B(:,9:11),2);
    
    AllSubjects_Data(iROI).Differential.MainEffect.DATA = C;
    AllSubjects_Data(iROI).Differential.MainEffect.MEAN = nanmean(C);
    AllSubjects_Data(iROI).Differential.MainEffect.STD = nanstd(C);
    AllSubjects_Data(iROI).Differential.MainEffect.SEM = nansem(C);
    
    clear B C
    
    % Compute main effects
    TEMP(:,1:3) = A(:,1:3);
    TEMP(:,end+1) = nanmean(A(:,1:3),2);
    TEMP(:,end+1:end+3) = A(:,4:6);
    TEMP(:,end+1) = nanmean(A(:,4:6),2);
    TEMP(:,end+1) = nanmean(A(:,[1 4]),2);
    TEMP(:,end+1) = nanmean(A(:,[2 5]),2);
    TEMP(:,end+1) = nanmean(A(:,[3 6]),2);
    TEMP(:,end+1) = nan(size(TEMP(:,end)));
    
    AllSubjects_Data(iROI).MainEffects.DATA = TEMP;
    AllSubjects_Data(iROI).MainEffects.MEAN = nanmean(TEMP);
    AllSubjects_Data(iROI).MainEffects.STD = nanstd(TEMP);
    AllSubjects_Data(iROI).MainEffects.SEM = nansem(TEMP);
    
    clear TEMP
    
    % Compute Bi - Uni
    TEMP = A(:,1:6);
    TEMP(:,end+1) = mean(cat(2,A(:,3,:)-A(:,1,:),A(:,6,:)-A(:,4,:)),2);
    TEMP(:,end+1) = mean(cat(2,A(:,3,:)-A(:,2,:),A(:,6,:)-A(:,5,:)),2);
    
    AllSubjects_Data(iROI).BiVSUni.DATA = TEMP;
    AllSubjects_Data(iROI).BiVSUni.MEAN = nanmean(TEMP);
    AllSubjects_Data(iROI).BiVSUni.STD = nanstd(TEMP);
    AllSubjects_Data(iROI).BiVSUni.SEM = nansem(TEMP);
    
    clear TEMP
    
     % Compute Interaction
     TEMP2(:,1,:) = A(:,3,:)-A(:,1,:) - (A(:,6,:)-A(:,4,:));
     TEMP2(:,2,:) = A(:,3,:)-A(:,2,:) - (A(:,6,:)-A(:,5,:));
     TEMP2(:,3,:) = A(:,2,:)-A(:,1,:) - (A(:,5,:)-A(:,4,:));
     
     AllSubjects_Data(iROI).Interaction.DATA = TEMP2;
     AllSubjects_Data(iROI).Interaction.MEAN = nanmean(TEMP2,3);
     AllSubjects_Data(iROI).Interaction.STD = nanstd(TEMP2,3);
     AllSubjects_Data(iROI).Interaction.SEM = nansem(TEMP2,3);
    
    
    
end




%% Saves
fprintf('\nSaving\n')

cd('/media/rxg243/BackUp2/AV_Integration_7T_2/Results/Profiles/Surfaces')
save(strcat('Data_BOLD_WholeROI.mat'))

cd(StartFolder)
