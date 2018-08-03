%%
% clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

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

Results_Folder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Results\Profiles\Surfaces';

CdtVec = repmat(1:6,12,1);
CdtVec = CdtVec(:)';

for iSub = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(iSub,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder,'BetaMapping','8Surf');
    
    
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
    fprintf(' Getting features A1:\n')
    
    iROI = 1;
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
    Cst_A1 = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_A1 = [Cst_A1;B(1,:)];
        end
    end
    
    Avg_A1 = squeeze(nanmean(nanmean(tmp,2),1));
    Avg_A1 = Avg_A1(:);

    [U_A1,~,~] = svd(Cst_A1);
    
    clear tmp Y X i B
    
    %%
    fprintf(' Getting features V1:\n')
    
    iROI = 3;
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
    Cst_V1 = [];
    Avg_V1 = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_V1 = [Cst_V1;B(1,:)];
            Avg_V1 = [Avg_V1; mean(Y)];
        end
    end
    
    [U_V1,~,~] = svd(Cst_V1);
    
    % split into 20 bins in case splitting by activated or deactivated
    % vertices is a bit too brutal
    V_logic = any([CdtVec==2 ; CdtVec==5]);
    [~,I_V1] = sort(nanmean(Cst_V1(V_logic,:)));
    IdxToAvg_V1 = floor(linspace(1,numel(I_V1),40+1));
    
    % split by activated and deactivated vertives
    Act_deact_V1 = nanmean(Cst_V1(V_logic,:))>0;
    
    clear tmp Y X i B
    
    clear Features Features_lh Features_rh
    
    % Stores Eigenvariates of A1 and V1
    Grp_U_V1(:,:,iSub) = U_V1;
    Grp_U_A1(:,:,iSub) = U_A1;
    
    
    
    %% compute connectivity
    fprintf(' Compute connectivity\n')
    
    for iSplit = 2:numel(IdxToAvg_V1)
        Vert_2_select = I_V1(IdxToAvg_V1(iSplit-1):IdxToAvg_V1(iSplit)); %specifies which vertex to take for that bin
        X = mean(mean(cat(3,Avg_V1(CdtVec==2,Vert_2_select),Avg_V1(CdtVec==5,Vert_2_select)),3),2);
        Y = mean([Avg_A1(CdtVec==2)  Avg_A1(CdtVec==5)],2);
        B_rank_avg(iSplit,iSub,:) = glmfit(X-nanmean(X), Y,'normal');
    end
    
    
    % CdtVec: column vector ((72 X 1) to identify which condition on which row
    % U_A1: output of SVD for A1 (72 by 72)
    % Cst_V1: Cst for each vertex of V1 for each block of each condition
    %   (72 X Nb vertices)
    % Act_deact_V1 : logical row vector (1 X Nb vertices) identifying is a vertex
    % is activated or deactivated by a V stim (averaged over attention
    % condition)
    
    % For each condition
    for iCdt = 1:6
        
        Row2Select = CdtVec==iCdt;
        
        Y = U_A1(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_A1(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        % On activated vertices
        X = Cst_V1(Row2Select,Act_deact_V1==1); % take the right condition S paramaters of V1
        X = nanmean(X,2); % Average across vertices (now a 1 X 12 vector)
        B_act(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
        B_act(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
        V1_Act(iSub) = nanmean(X); %Keep the mean activation across blocks
        
        % On deactivated vertices same as above
        X = Cst_V1(Row2Select,Act_deact_V1==0);
        X = nanmean(X,2);
        B_deact(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal');
        B_deact(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal');
        V1_Deact(iSub) = nanmean(X);
        
        % Same as above but for different groups of vertices with
        % increasing level of activation in V1
        % I_V1: defines which vertices to take for each "bin"
        for iSplit = 2:numel(IdxToAvg_V1)
            Vert_2_select = I_V1(IdxToAvg_V1(iSplit-1):IdxToAvg_V1(iSplit)); %specifies which vertex to take for that bin
            X = Cst_V1(Row2Select,Vert_2_select);
            X = nanmean(X,2);
            B_rank(iSplit,iSub,iCdt,:,1) = glmfit(X-nanmean(X), Y,'normal');
            B_rank(iSplit,iSub,iCdt,:,2) = glmfit(X-nanmean(X), Y2,'normal');
            V1_Act_rank(iSub,iSplit,iCdt) = nanmean(X);
        end
        
    end
    
    
    % Pooled over A and V only separately for each attention condition
    for iCdt = 1:2
        
        if iCdt==1
            Row2Select = ismember(CdtVec,1:2);
        else
            Row2Select = ismember(CdtVec,4:5);
        end
        
        Y = U_A1(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_A1(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        % On activated vertices
        X = Cst_V1(Row2Select,Act_deact_V1==1); % take the right condition S paramaters of V1
        X = nanmean(X,2); % Average across vertices (now a 1 X 12 vector)
        B_act_att(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
        B_act_att(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
        V1_Act_att(iSub) = nanmean(X); %Keep the mean activation across blocks
        
        % On deactivated vertices same as above
        X = Cst_V1(Row2Select,Act_deact_V1==0);
        X = nanmean(X,2);
        B_deact_att(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal');
        B_deact_att(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal');
        V1_Deact_att(iSub) = nanmean(X);
        
        % Same as above but for different groups of vertices with
        % increasing level of activation in V1
        % I_V1: defines which vertices to take for each "bin"
        for iSplit = 2:numel(IdxToAvg_V1)
            Vert_2_select = I_V1(IdxToAvg_V1(iSplit-1):IdxToAvg_V1(iSplit)); %specifies which vertex to take for that bin
            X = Cst_V1(Row2Select,Vert_2_select);
            X = nanmean(X,2);
            B_rank_att(iSplit,iSub,iCdt,:,1) = glmfit(X-nanmean(X), Y,'normal');
            B_rank_att(iSplit,iSub,iCdt,:,2) = glmfit(X-nanmean(X), Y2,'normal');
            V1_Act_rank_att(iSub,iSplit) = nanmean(X);
        end
        
    end
    
    
    % Pooled over A att and V att separately for each V and A condition
    for iCdt = 1:2
        
        if iCdt==1
            Row2Select = ismember(CdtVec,[1 4]);
        else
            Row2Select = ismember(CdtVec,[2 5]);
        end
        
        Y = U_A1(Row2Select,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
        Y2 = nanmean(Cst_A1(Row2Select,:),2); % Average across vertices of activation (now a 1 X 12 vector)
        
        % On activated vertices
        X = Cst_V1(Row2Select,Act_deact_V1==1); % take the right condition S paramaters of V1
        X = nanmean(X,2); % Average across vertices (now a 1 X 12 vector)
        B_act_sens(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
        B_act_sens(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
        V1_Act_sens(iSub) = nanmean(X); %Keep the mean activation across blocks
        
        % On deactivated vertices same as above
        X = Cst_V1(Row2Select,Act_deact_V1==0);
        X = nanmean(X,2);
        B_deact_sens(iSub,iCdt,:,1) = glmfit(X-nanmean(X),Y,'normal');
        B_deact_sens(iSub,iCdt,:,2) = glmfit(X-nanmean(X),Y2,'normal');
        V1_Deact_sens(iSub) = nanmean(X);
        
        % Same as above but for different groups of vertices with
        % increasing level of activation in V1
        % I_V1: defines which vertices to take for each "bin"
        for iSplit = 2:numel(IdxToAvg_V1)
            Vert_2_select = I_V1(IdxToAvg_V1(iSplit-1):IdxToAvg_V1(iSplit)); %specifies which vertex to take for that bin
            X = Cst_V1(Row2Select,Vert_2_select);
            X = nanmean(X,2);
            B_rank_sens(iSplit,iSub,iCdt,:,1) = glmfit(X-nanmean(X), Y,'normal');
            B_rank_sens(iSplit,iSub,iCdt,:,2) = glmfit(X-nanmean(X), Y2,'normal');
            V1_Act_rank_sens(iSub,iSplit) = nanmean(X);
        end
        
    end
    
    
    cd(Results_Folder)
    
    clear Cst_V1 Cst_A1 iCdt iSplit
    
    
end

%%

save(fullfile(Results_Folder,strcat('Connectivity_Surf_', num2str(NbLayers), '_layers.mat')), ...
    'Grp_U_V1', 'Grp_U_A1', 'B_rank_avg',...
    'B_rank', 'B_act', 'B_deact', ...
    'V1_Act', 'V1_Deact', 'V1_Act_rank',...
    'B_rank_att', 'B_act_att', 'B_deact_att', ...
    'V1_Act_att', 'V1_Deact_att', 'V1_Act_rank_att',...
    'B_rank_sens', 'B_act_sens', 'B_deact_sens', ...
    'V1_Act_sens', 'V1_Deact_sens', 'V1_Act_rank_sens')


