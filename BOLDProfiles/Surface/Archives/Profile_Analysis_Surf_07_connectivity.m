%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

Print=0;

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

SubjectList = [...
    %         '02';...
    %         '03';...
    '04';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    %         '13';...
    '15';...
    '16'
    ];


MinLayer = 1;
NbLayers = 7;
NbLayers = NbLayers+1;
Ind = NbLayers:-1:1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
% DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);

Results_Folder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Results\Profiles\Surfaces';


for SubjInd = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
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
    
    % Format for reading the vertices from the VTK file
    Spec = repmat('%f ', 1, NbLayers);
    
    NbVertices = nan(1,2);
    
    % For the 2 hemispheres
    clear VertexWithDataHS MappingBothHS
    
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
        
        %                 InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
        %                     ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
        %
        %                 [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);
        %                 NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData;
            MappingBothHS{hs} = AllMapping;
        else
            
            AllMapping = nan(NbVertices(hs),NbLayers,size(Betas,1));
            
            fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
            
            for iBeta = 1:size(Betas,1)
                
                A = fileread(fullfile(Data_Folder,Betas(iBeta).name)); % reads file quickly
                B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
                
                C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
                Mapping = cell2mat(C); %clear C
                
                if size(Mapping,1)~=(NbVertices(hs)) %#ok<PFBNS>
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
            
            
            save(FeatureSaveFile,'Vertex','Face','AllMapping','VertexWithData', '-v7.3')
            
            VertexWithDataHS{hs} = VertexWithData;
            MappingBothHS{hs} = AllMapping;
        end
        
        clear Betas InfSurfFile Mapping A B C
        
    end
    
    %         if any(NbVertex ~= NbVertices)
    %             NbVertex
    %             NbVertices
    %             error('The number of vertices does not match.')
    %         end
    
    Features_lh = nan(NbVertex(1),NbLayers,size(MappingBothHS{1},3));
    Features_lh(VertexWithDataHS{1},:,:) = MappingBothHS{1};
    
    Features_rh = nan(NbVertex(2),NbLayers,size(MappingBothHS{2},3));
    Features_rh(VertexWithDataHS{2},:,:) = MappingBothHS{2};
    
    %%
    fprintf(' Getting features V1:\n')
    
    iROI = 3;
    Data_ROI.name = ROI(iROI).name;
    fprintf(['  '  Data_ROI.name '\n'])
    
    Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
        Features_rh(ROI(iROI).VertOfInt{2},:,:));
    
    % get only V activation
    Cdts = [2 5];
    Beta2Sel = [];
    for BlockInd = 1:3
        Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
            [Conditions_Names{Cdts(1)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
    end
    Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
    tmp = Features(:,2:7,Beta2Sel);
    
    Beta2Sel = [];
    for BlockInd = 1:3
        Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
            [Conditions_Names{Cdts(2)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
    end
    Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
    tmp(:,:,:,2) = Features(:,2:7,Beta2Sel);
    
    tmp = mean(tmp,4); %averages across attention conditions
    tmp(any(any(isnan(tmp),2),3),:,:) = []; % remove nans
    tmp(any(any(tmp==0,2),3),:,:) = []; % remove 0 values from masking
    
    % run a laminar GLM to get the constant to sort vertices
    X=repmat(DesMat,size(tmp,3),1);
    Y = shiftdim(tmp,1);
    Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
    Cst_lin = pinv(X)*Y;
    
    % split into 20 bins in case splitting by activated or deactivated
    % vertices is a bit too brutal
    [~,I] = sort(Cst_lin(1,:));
    IdxToAvg = floor(linspace(1,numel(I),20+1));
    for iSplit = 2:numel(IdxToAvg)
        Vert_2_select = I(IdxToAvg(iSplit-1):IdxToAvg(iSplit));
        V1_act_rank_V_layers(:,:,iSplit-1) = squeeze(mean(tmp(Vert_2_select,:,:)));
        V1_act_rank_V(iSplit-1,:) = mean(V1_act_rank_V_layers(:,:,iSplit-1)); % average across layers
        V1_act_rank_beta(iSplit-1,SubjInd) = mean(V1_act_rank_V(iSplit-1,:));
    end
    
    % split by activated and deactivated vertives
    Act_deact = Cst_lin(1,:)>0;
    V1_Act_V_layers = squeeze(mean(tmp(Act_deact,:,:)));
    V1_Deact_V_layers = squeeze(mean(tmp(Act_deact==0,:,:)));
    
     
    
    % averages across layers
    tmp = squeeze(mean(tmp,2));
    V1_V_Stim_Vert = tmp';
    V1_Act_V = mean(tmp(Act_deact,:));
    V1_Deact_V = mean(tmp(Act_deact==0,:));
    
    clear tmp
    
    %%
    fprintf(' Getting features A1:\n')
    
    iROI = 1;
    Data_ROI.name = ROI(iROI).name;
    fprintf(['  '  Data_ROI.name '\n'])
    
    Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
        Features_rh(ROI(iROI).VertOfInt{2},:,:));
    
    % get only V deactivation
    Cdts = [2 5];
    Beta2Sel = [];
    for BlockInd = 1:3
        Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
            [Conditions_Names{Cdts(1)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
    end
    Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
    tmp = Features(:,2:7,Beta2Sel);
    
    Beta2Sel = [];
    for BlockInd = 1:3
        Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
            [Conditions_Names{Cdts(2)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
    end
    Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
    tmp(:,:,:,2) = Features(:,2:7,Beta2Sel);
    
    tmp = mean(tmp,4); %averages across attention conditions
    tmp(any(any(isnan(tmp),2),3),:,:) = []; % remove nans
    tmp(any(any(tmp==0,2),3),:,:) = []; % remove 0 values from masking
    
    % average across vertices
    A1_Deact_V_layers = squeeze(mean(tmp));
    
    % averages across layers
    tmp = squeeze(mean(tmp,2));
    A1_V_Stim_Vert = tmp';
    A1_Deact_V = mean(tmp);
    
    clear tmp;
    
    %% compute connectivity
    % from layer to layer
    for iLayer = 1:6
        for jLayer = 1:6
            % for activated vertices of V1
            B_act_layers(iLayer,jLayer,:,SubjInd) = glmfit(V1_Act_V_layers(iLayer,:)-mean(V1_Act_V_layers(iLayer,:)), ...
                A1_Deact_V_layers(jLayer,:),'normal');
            
            % for activated vertices of V1
            B_deact_layers(iLayer,jLayer,:,SubjInd) = glmfit(V1_Deact_V_layers(iLayer,:)-mean(V1_Deact_V_layers(iLayer,:)), ...
                A1_Deact_V_layers(jLayer,:),'normal');
            
            % with 20 bins of activation level
            for iSplit = 1:20
                B_rank_layers(iLayer,jLayer,iSplit,:,SubjInd) = glmfit(V1_act_rank_V_layers(iLayer,:,iSplit)-mean(V1_act_rank_V_layers(iLayer,:,iSplit)), ...
                    A1_Deact_V_layers(jLayer,:),'normal');
            end
        end
    end
    
    % same but without the layer aspect
    B_act(SubjInd,:) = glmfit(V1_Act_V-mean(V1_Act_V), A1_Deact_V ,'normal');
    B_deact(SubjInd,:) = glmfit(V1_Deact_V-mean(V1_Deact_V), A1_Deact_V,'normal');
    for iSplit = 1:20
        B_rank(iSplit,:,SubjInd) = glmfit(V1_act_rank_V(iSplit,:)-mean(V1_act_rank_V(iSplit,:)), A1_Deact_V,'normal');
    end
    
    
    %% Compute SVD
%     [U, S, V] = svd(A1_V_Stim_Vert'*V1_V_Stim_Vert);
end

%%
cd(Results_Folder)
save(strcat('Connectivity_Surf_Block_V1_to_A1_', num2str(NbLayers), '_layers_old.mat'), ...
    'B_rank', 'B_act', 'B_deact', 'B_rank_layers', 'B_act_layers', 'B_deact_layers', 'V1_act_rank_beta')

return

%%
close all
FigDim = [100 100 1000 700];
figure('name', 'rank', 'Position', FigDim, 'Color', [1 1 1]);

% subplot(121)
% hold on
% grid on
% herrorbar(mean(B_rank(:,1,:),3),1:20,nansem(B_rank(:,1,:),3),'o--k')
% t=herrorbar(mean(B_act(:,1)),23,nansem(B_act(:,1)),'.k')
% set(t, 'linewidth', 2, 'markersize', 20)
% t=herrorbar(mean(B_deact(:,1)),22,nansem(B_act(:,1)),'.k')
% set(t, 'linewidth', 2, 'markersize', 20)
% 
% ylabel('[V-Fix]')
% set(gca,'ygrid','off','xtick',-3:.25:3,'xticklabel',-3:.25:3, ...
%     'ytick', [1 20 22 23],  'yticklabel',{'more deactivated';'more activated'; 'deactivated'; 'activated'})
% title('Connectivity ; intercept')
% axis([-1.25 1 0.5 23.5])
% 
% subplot(122)
hold on
grid on

% t=herrorbar(mean(B_rank(:,2,:),3),mean(V1_act_rank_beta,2),nansem(B_rank(:,2,:),3),'o-b')
% set(t, 'linewidth', 1.5, 'markersize', 5)
for iSub = 1:size(B_rank,3)
    plot(B_rank(:,2,iSub),V1_act_rank_beta(:,iSub),'color','k') % [.3 .3 1]
end

t=errorbar(mean(B_rank(:,2,:),3),mean(V1_act_rank_beta,2),....
    nansem(V1_act_rank_beta,2),nansem(V1_act_rank_beta,2),...
    nansem(B_rank(:,2,:),3),nansem(B_rank(:,2,:),3),'o-b')
set(t, 'linewidth', 1.5, 'markersize', 5)

t=herrorbar(mean(B_act(:,2)),15,nansem(B_act(:,2)),'.b')
set(t, 'linewidth', 3, 'markersize', 25)
t=herrorbar(mean(B_deact(:,2)),14,nansem(B_act(:,2)),'.b')
set(t, 'linewidth', 3, 'markersize', 25)
plot([-50 50],[0 0], '--k')
plot([0 0],[-50 50], '--k')
ylabel('[V-Fix]')
set(gca,'ygrid','off','xtick',-3:.1:3,'xticklabel',-3:.1:3, ...
    'ytick', [-4:2:12 14 15],  'yticklabel',{'-4';'-2';'0';'2';'4';'6';'8';'10';'12'; 'deactivated'; 'activated'})
title('Connectivity ; slope')
axis([-.5 1 -6 16])
% axis([0 .4 -4 15])

print(gcf, fullfile(FigureFolder,'Baseline', 'Connectvity_V1_to_A1_rank.tif'), '-dtiff')


