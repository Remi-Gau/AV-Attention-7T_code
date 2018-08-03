%%
clc; clear;

V1_2_A1 = 0;
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

Session_Vec = repmat((1:4)',3*6,1);

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
    
    Avg_TargetROI = squeeze(nanmean(nanmean(tmp,2),1));
    Avg_TargetROI = Avg_TargetROI(:);
    
    
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
    Avg_SeedROI = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_SeedROI = [Cst_SeedROI;B(1,:)];
            Avg_SeedROI = [Avg_SeedROI; mean(Y)];
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
    
    % split into 20 bins in case splitting by activated or deactivated
    % vertices is a bit too brutal
    if V1_2_A1==0
        Logic = any([CdtVec==1 ; CdtVec==4]);
    else
        Logic = any([CdtVec==2 ; CdtVec==5]);
    end

    clear I_V1
    [~,I_V1(1,:)] = sort(nanmean(Cst_SeedROI(all([Logic' ismember(Session_Vec,1:2)],2),:)));
    [~,I_V1(2,:)] = sort(nanmean(Cst_SeedROI(all([Logic' ismember(Session_Vec,3:4)],2),:)));
    IdxToAvg = floor(linspace(1,size(I_V1,2),40+1));
    
    % split by activated and deactivated vertives
    clear Act_deact_Seed
    Act_deact_Seed(1,:) = nanmean(Cst_SeedROI(all([Logic' ismember(Session_Vec,1:2)],2),:))>0;
    Act_deact_Seed(2,:) = nanmean(Cst_SeedROI(all([Logic' ismember(Session_Vec,3:4)],2),:))>0;
    
    clear tmp Y X i B
    
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
    for iCV = 1:2
        if iCV==1
            test_ses = ismember(Session_Vec,3:4);
        else
            test_ses = ismember(Session_Vec,1:2);
        end

        % For each condition
        for iCdt = 1:6
            
            Row2Select = (CdtVec==iCdt)';

            Test_rows = all([Row2Select test_ses],2);
            
            Y = U_TargetROI(Test_rows,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
            Y2 = nanmean(Cst_TargetROI(Test_rows,:),2); % Average across vertices of activation (now a 1 X 12 vector)
            
            % Run GLM on all vertices
            for iVert = 1:size(Cst_SeedROI,2)
                X = Cst_SeedROI(Test_rows,iVert); % take the right condition S paramaters of V1
                Vert_BOLD(:,1,iVert) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
                Vert_BOLD(:,2,iVert) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
            end
            
            % On activated vertices
            B_act(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==1),3);
            B_act(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==1),3);
            Seed_Act(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==1))); %Keep the mean activation across blocks
            
            % On deactivated vertices same as above
            B_deact(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==0),3);
            B_deact(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==0),3);
            Seed_Deact(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==0)));
            
            % Same as above but for different groups of vertices with
            % increasing level of activation in V1
            % I_V1: defines which vertices to take for each "bin"
            for iSplit = 2:numel(IdxToAvg)
                Vert_2_select = I_V1(iCV,IdxToAvg(iSplit-1):IdxToAvg(iSplit)); %specifies which vertex to take for that bin
                B_rank(iSplit,iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Vert_2_select),3);
                B_rank(iSplit,iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Vert_2_select),3);
                Seed_Act_rank(iSub,iSplit,iCdt,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Vert_2_select)));
            end
            
            % Taking all vertices
            B_all(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,:),3);
            B_all(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,:),3);
            
            
        end
        
        
        % Pooled over A and V only separately for each attention condition
        for iCdt = 1:2
            
            if iCdt==1
                Row2Select = ismember(CdtVec,1:2)';
            else
                Row2Select = ismember(CdtVec,4:5)';
            end
            
            Test_rows = all([Row2Select test_ses],2);
            
            Y = U_TargetROI(Test_rows,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
            Y2 = nanmean(Cst_TargetROI(Test_rows,:),2); % Average across vertices of activation (now a 1 X 12 vector)
            
            for iVert = 1:size(Cst_SeedROI,2)
                X = Cst_SeedROI(Test_rows,iVert); % take the right condition S paramaters of V1
                Vert_BOLD(:,1,iVert) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
                Vert_BOLD(:,2,iVert) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
            end
            
            % On activated vertices
            B_act_att(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==1),3);
            B_act_att(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==1),3);
            Seed_Act_att(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==1)));
            
            % On deactivated vertices same as above
            B_deact_att(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==0),3);
            B_deact_att(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==0),3);
            Seed_Deact_att(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==0)));
            
            % Same as above but for different groups of vertices with
            % increasing level of activation in V1
            % I_V1: defines which vertices to take for each "bin"
            for iSplit = 2:numel(IdxToAvg)
                Vert_2_select = I_V1(iCV,IdxToAvg(iSplit-1):IdxToAvg(iSplit)); %specifies which vertex to take for that bin
                B_rank_att(iSplit,iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Vert_2_select),3);
                B_rank_att(iSplit,iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Vert_2_select),3);
                Seed_Act_rank_att(iSub,iSplit,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Vert_2_select)));
            end
            
            % Taking all vertices
            B_all_att(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,:),3);
            B_all_att(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,:),3);
            
        end
        
        
        % Pooled over A att and V att separately for each V and A condition
        for iCdt = 1:2
            
            if iCdt==1
                Row2Select = ismember(CdtVec,[1 4])';
            else
                Row2Select = ismember(CdtVec,[2 5])';
            end

            Test_rows = all([Row2Select test_ses],2);
            
            Y = U_TargetROI(Test_rows,1); % take the correct condition of the first Eigenvariate of A1 (now a 1 X 12 vector)
            Y2 = nanmean(Cst_TargetROI(Test_rows,:),2); % Average across vertices of activation (now a 1 X 12 vector)
            
            for iVert = 1:size(Cst_SeedROI,2)
                X = Cst_SeedROI(Test_rows,iVert); % take the right condition S paramaters of V1
                Vert_BOLD(:,1,iVert) = glmfit(X-nanmean(X),Y,'normal'); % Regress (with mean centering of V1 values) on Eigenvariate
                Vert_BOLD(:,2,iVert) = glmfit(X-nanmean(X),Y2,'normal'); % Same but on activation levels
            end
            
            % On activated vertices
            B_act_sens(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==1),3);
            B_act_sens(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==1),3);
            Seed_Act_sens(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==1)));
            
            % On deactivated vertices same as above
            B_deact_sens(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Act_deact_Seed(iCV,:)==0),3);
            B_deact_sens(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Act_deact_Seed(iCV,:)==0),3);
            Seed_Deact_sens(iSub,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Act_deact_Seed(iCV,:)==0)));
            
            % Same as above but for different groups of vertices with
            % increasing level of activation in V1
            % I_V1: defines which vertices to take for each "bin"
            for iSplit = 2:numel(IdxToAvg)
                Vert_2_select = I_V1(iCV,IdxToAvg(iSplit-1):IdxToAvg(iSplit)); %specifies which vertex to take for that bin
                B_rank_sens(iSplit,iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,Vert_2_select),3);
                B_rank_sens(iSplit,iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,Vert_2_select),3);
                Seed_Act_rank_sens(iSub,iSplit,iCV) = nanmean(nanmean(Cst_SeedROI(Test_rows,Vert_2_select)));
            end
            
            % Taking all vertices
            B_all_sens(iSub,iCdt,:,1,iCV) = nanmean(Vert_BOLD(:,1,:),3);
            B_all_sens(iSub,iCdt,:,2,iCV) = nanmean(Vert_BOLD(:,2,:),3);
            
        end
        
    end
    
    cd(Results_Folder)
    
    clear Cst_SeedROI Cst_TargetROI iCdt iSplit
    
    
end

%%

save(fullfile(Results_Folder,strcat('Connectivity_', Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')), ...
    'Grp_U_Target', 'Grp_U_Seed',...
    'B_rank', 'B_act', 'B_deact', ...
    'B_all', 'B_all_att', 'B_all_sens', ...
    'Seed_Act', 'Seed_Deact', 'Seed_Act_rank',...
    'B_rank_att', 'B_act_att', 'B_deact_att', ...
    'Seed_Act_att', 'Seed_Deact_att', 'Seed_Act_rank_att',...
    'B_rank_sens', 'B_act_sens', 'B_deact_sens', ...
    'Seed_Act_sens', 'Seed_Deact_sens', 'Seed_Act_rank_sens')