%%
clc; clear;

V1_2_A1 = 1;
if V1_2_A1==1
    Suffix = 'V1_2_A1';
else
    Suffix = 'A1_2_V1';
end

addpath(genpath(fullfile('/home','rxg243','Dropbox',...
    'Code','MATLAB','Neuroimaging','SPM','spm12')))

StartDirectory = fullfile(pwd, '..', '..', '..');
cd(StartDirectory)
addpath(genpath(fullfile(StartDirectory, 'SubFun')))
Get_dependencies('/home/rxg243/Dropbox')

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

MinLayer = 1;
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
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);

Results_Folder = fullfile(StartDirectory,'Results');
mkdir(Results_Folder)

CdtVec = repmat(1:6,12,1);
CdtVec = CdtVec(:)';

for iSub = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(iSub,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartDirectory, ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder);
    
    Data_Folder = fullfile(SubjectFolder);
    
    
    %% Creates a cell that lists the names of the beta images as well as their column number in the design matrix
    load(fullfile(SubjectFolder, 'SPM.mat'))
    
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
    load(fullfile(SubjectFolder,['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
    %% Read features
    fprintf(' Reading VTKs\n')
    
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
        
        
        %% Load data or extract them
        load(FeatureSaveFile)
        VertexWithDataHS{hs} = VertexWithData;
        MappingBothHS{hs} = AllMapping;
        
    end
    
    Features_lh = nan(NbVertex(1),NbLayers,size(MappingBothHS{1},3));
    Features_lh(VertexWithDataHS{1},:,:) = MappingBothHS{1};
    
    Features_rh = nan(NbVertex(2),NbLayers,size(MappingBothHS{2},3));
    Features_rh(VertexWithDataHS{2},:,:) = MappingBothHS{2};
    
    clear MappingBothHS AllMapping VertexWithData Face Vertex VertexWithDataHS NbVertex
    
    %%
    fprintf(' Getting features Target:\n')
    
    if V1_2_A1==0
        iROI = 3;
    else
        iROI = 1;
    end
    Data_ROI.name = ROI(iROI).name;
    fprintf(['  '  Data_ROI.name '\n'])
    
    Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
        Features_rh(ROI(iROI).VertOfInt{2},:,:));
    
    % get all conditions
    tmp = [];
    for iCdt = 1:6
        Beta2Sel = [];
        for BlockInd = 1:3
            Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                [Conditions_Names{iCdt} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
        end
        Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
        tmp(:,:,:,iCdt) = Features(:,2:7,Beta2Sel);
    end
    tmp(any(any(any(isnan(tmp),2),3),4),:,:,:) = []; % remove nans
    tmp(any(any(any(tmp==0,2),3),4),:,:,:) = []; % remove 0 values from masking
    
    % run a laminar GLM to get the constant
    Cst_TargetROI = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_TargetROI = [Cst_TargetROI;B(1,:)];
        end
    end
    
    % copied from spm_regions to ensure polarity of first eigen variate
    [m,n]   = size(Cst_TargetROI);
    [u,s,u] = svd(Cst_TargetROI*Cst_TargetROI');
    s       = diag(s);
    u       = u(:,1);
    v       = Cst_TargetROI'*u/sqrt(s(1));
    
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    U_TargetROI       = u*sqrt(s(1)/n);
    clear s u v d
    
    clear tmp Y X i B
    
    %%
    fprintf(' Getting features Seed:\n')
    
    if V1_2_A1==0
        iROI = 1;
    else
        iROI = 3;
    end
    Data_ROI.name = ROI(iROI).name;
    fprintf(['  '  Data_ROI.name '\n'])
    
    Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
        Features_rh(ROI(iROI).VertOfInt{2},:,:));
    
    % get all conditions
    tmp = [];
    for iCdt = 1:6
        Beta2Sel = [];
        for BlockInd = 1:3
            Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                [Conditions_Names{iCdt} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
        end
        Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
        tmp(:,:,:,iCdt) = Features(:,2:7,Beta2Sel);
    end
    tmp(any(any(any(isnan(tmp),2),3),4),:,:,:) = []; % remove nans
    tmp(any(any(any(tmp==0,2),3),4),:,:,:) = []; % remove 0 values from masking
    
    % run a laminar GLM to get the constant
    Cst_SeedROI = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_SeedROI = [Cst_SeedROI;B(1,:)];
        end
    end
    
    % copied from spm_regions to ensure polarity of first eigen variate
    [m,n]   = size(Cst_SeedROI);
    [u,s,u] = svd(Cst_SeedROI*Cst_SeedROI');
    s       = diag(s);
    u       = u(:,1);
    v       = Cst_SeedROI'*u/sqrt(s(1));
    
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    U_SeedROI       = u*sqrt(s(1)/n);
    clear s u v d
    
    clear Features Features_lh Features_rh
    
    % Stores Eigenvariates of A1 and V1
    Grp_U_Seed(:,:,iSub) = U_SeedROI;
    Grp_U_Target(:,:,iSub) = U_TargetROI;
    
    %% compute connectivity
    fprintf(' Compute connectivity\n')
    % CdtVec: column vector ((72 X 1) to identify which condition on which row
    % U_A1: output of SVD for A1 (72 by 72)
    % Cst_V1: Cst for each vertex of V1 for each block of each condition
    %   (72 X Nb vertices)
    % Act_deact_V1 : logical row vector (1 X Nb vertices) identifying is a vertex
    % is activated or deactivated by a V stim (averaged over attention
    % condition)
    Grp_Cst_Seed(:,iSub) = nanmean(Cst_TargetROI,2);
    Grp_Cst_Target(:,iSub) = nanmean(Cst_SeedROI,2);
    
    % For each condition
    for iCdt = 1:6
        
        Row2Select = (CdtVec==iCdt)';
        
        Y = U_TargetROI(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_TargetROI(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        X = nanmean(Cst_SeedROI(Row2Select,:),2); % take the right condition S paramaters of V1
%         B_all(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
%         B_all(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels

        B_all(iSub,iCdt,:,1) = regress(Y,X); % Regress (with mean centering of V1 values) on Eigenvariate
        B_all(iSub,iCdt,:,2) = regress(Y2,X); % Same but on activation levels
        
    end
    
    
    % Pooled over A and V only separately for each attention condition
    for iCdt = 1:2
        
        if iCdt==1
            Row2Select = ismember(CdtVec,1:2)';
        else
            Row2Select = ismember(CdtVec,4:5)';
        end
        
        Y = U_TargetROI(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_TargetROI(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        X = nanmean(Cst_SeedROI(Row2Select,:),2); % take the right condition S paramaters of V1
%         B_all_att(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
%         B_all_att(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
        
        B_all_att(iSub,iCdt,:,1) = regress(Y,X); % Regress (with mean centering of V1 values) on Eigenvariate
        B_all_att(iSub,iCdt,:,2) = regress(Y2,X); % Same but on activation levels
        
        
    end
    
    
    % Pooled over A att and V att separately for each V and A condition
    for iCdt = 1:2
        
        if iCdt==1
            Row2Select = ismember(CdtVec,[1 4])';
        else
            Row2Select = ismember(CdtVec,[2 5])';
        end
        
        Y = U_TargetROI(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_TargetROI(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        X = nanmean(Cst_SeedROI(Row2Select,:),2); % take the right condition S paramaters of V1
%         B_all_sens(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
%         B_all_sens(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
        
        B_all_sens(iSub,iCdt,:,1) = regress(Y,X); % Regress (with mean centering of V1 values) on Eigenvariate
        B_all_sens(iSub,iCdt,:,2) = regress(Y2,X); % Same but on activation levels
  
    end

    cd(Results_Folder)
    
    clear Cst_SeedROI Cst_TargetROI iCdt iSplit
    
    
end

%%

save(fullfile(Results_Folder,strcat('Connectivity_all_', Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')), ...
    'B_all', 'B_all_att', 'B_all_sens',...
    'Grp_Cst_Seed', 'Grp_Cst_Target', 'Grp_U_Seed', 'Grp_U_Target')