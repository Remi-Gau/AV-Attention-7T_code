function Profile_Analysis_Surf_05_Cluster
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

NbLayers = 7;
NbLayers = NbLayers+1;





NbWorkers = 4;

MatlabVer = version('-release');
if str2double(MatlabVer(1:4))>2013
    if isempty(gcp)
        KillGcpOnExit = 1;
        parpool(NbWorkers);
    else
        KillGcpOnExit = 0;
    end
    % else
    %     if matlabpool('size') == 0
    %         KillGcpOnExit = 1;
    %         matlabpool(NbWorkers)
    %     elseif matlabpool('size') ~= NbWorkers
    %         matlabpool close
    %         matlabpool(NbWorkers)
    %         KillGcpOnExit = 0;
    %     else
    %         KillGcpOnExit = 0;
    %     end
end





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
        
        close all
        
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
        
        [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);
        NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData;
        else
            
            AllMapping = nan(NbVertices(hs),NbLayers,size(Betas,1));
            
            fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
            
            parfor iBeta = 1:size(Betas,1)
                
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
            
            VertexWithDataHS{hs} = VertexWithData;
            
            save(FeatureSaveFile,'Vertex','Face','AllMapping','VertexWithData', '-v7.3')
        end
        
        clear Betas InfSurfFile Mapping A B C
        
        
        %% Run for A, V and AV
        for CondInd = 1:length(Conditions_Names)/2 % For each Condition
            
            Beta2Sel = [];
            Beta2Sel2 = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
                Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd+3} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
            
            Features = mean(cat(4,AllMapping(:,2:end-1,Beta2Sel), AllMapping(:,2:end-1,Beta2Sel2)),4);
            if sum(isnan(Features(:)))>0
                warning('We have %i NaNs for %s', sum(isnan(Features(:))), Conditions_Names{CondInd})
            end
            if sum(Features(:)==0)>0
                warning('We have %i zeros for %s', sum(Features(:)==0))
            end
            
            Features = mean(Features,3);
            
            X = nan(NbVertex(hs),NbLayers-2);
            X(VertexWithDataHS{hs},:) = Features;
            
            Name{1} = fullfile(Results_Folder,'Cdtions',['Subj_' SubjID '_' ]);
            Name{2} = ['_' HsSufix 'cr_' strrep(Conditions_Names{CondInd},' ','')];
        
            Clustering_BOLD(ROI, X , hs, NbVertex(hs), Vertex, Face, Name)
            
            clear Features Beta2Sel X Name
            
        end
        
        
        %% Run GLMs for Attention
        Cond2Contrast = {...
            1, 4;...
            2, 5;...
            3, 6;...
            };
        
        AttCondNames = {...
            'A_Stim_A_Att-V_Att'; ...
            'V_Stim_A_Att-V_Att'; ...
            'AV_Stim_A_Att-V_Att'};
        
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
            
            FeaturesAll(:,:,:,AttCondInd) = Features; %#ok<*SAGROW>
            
            X = nan(NbVertex(hs),NbLayers-2);
            X(VertexWithDataHS{hs},:) = Features;

            Name{1} = fullfile(Results_Folder,'Att',['Subj_' SubjID '_' ]);
            Name{2} = ['_' HsSufix 'cr_' strrep(AttCondNames{AttCondInd},' ','')];
            
            Clustering_BOLD(ROI, X , hs, NbVertex(hs), Vertex, Face, Name)
            
            clear Features Beta2Sel X Name
            
        end
        
        Features = mean(FeaturesAll,4);
        
        X = nan(NbVertex(hs),NbLayers-2);
        X(VertexWithDataHS{hs},:) = Features;

        Name{1} = fullfile(Results_Folder,'Att',['Subj_' SubjID '_' ]);
        Name{2} = ['_' HsSufix 'cr_A_Att-V_Att'];
        
        Clustering_BOLD(ROI, X , hs, NbVertex(hs), Vertex, Face, Name)
        
        clear Features X Name
        
               
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
                        [Conditions_Names{Cond2Contrast{CrossSensCondInd,1}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
                    Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                        [Conditions_Names{Cond2Contrast{CrossSensCondInd,2}(CondInd)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>.
                end
                Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
                Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
                
                Features(:,:,:,CondInd) = AllMapping(:,2:end-1,Beta2Sel) - ...
                    AllMapping(:,2:end-1,Beta2Sel2);
            end
            
            Features = mean(Features,4);
            
            X = nan(NbVertex(hs),NbLayers-2);
            X(VertexWithDataHS{hs},:) = Features;

            Name{1} = fullfile(Results_Folder,'CrossSens',['Subj_' SubjID '_' ]);
            Name{2} = ['_' HsSufix 'cr_' strrep(CrossSensCondNames{CrossSensCondInd},' ','')];
            
            Clustering_BOLD(ROI, X , hs, NbVertex(hs), Vertex, Face, Name)
            
            clear Features Beta2Sel X Name
            
        end

        
    end
    
    
    
    clear Cdt Name ToPlot Range Cdt Name ToPlot Range BetaCdtion BetaCrossSens BetaAtt
    
    
end

if str2double(MatlabVer(1:4))>2013
    if KillGcpOnExit
        delete(gcp)
    end
else
end

% save(fullfile(StartFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'), ...
%     'ROI', 'GrpBetaStimSpeAtt', 'GrpBetaStimAtt', 'GrpCrossModAtt', 'GrpBetaCrossSens','GrpBetaAtt','GrpBetaCdt')

cd(StartFolder)

end



function Clustering_BOLD(ROI, X , hs, NbVertex, Vertex, Face, Name)

MaxNbCluster = 30;
Replicates = 100;
MaxIter = 200;

opts = statset('MaxIter', MaxIter, 'UseParallel', true);


FigDim = [100, 100, 1000, 1500];
FontSize = 10;

Visibility = 'on';

ClusterColors = [...
    166,206,227;...
    0,0,0;...
    31,120,180;...
    178,223,138;...
    51,160,44;...
    251,154,153;...
    227,26,28;...
    253,191,111;...
    255,127,0;...
    202,178,214;...
    106,61,154;...
    255,255,153;...
    177,89,40;...
    228,26,28;...
    55,126,184;...
    77,175,74;...
    152,78,163;...
    255,127,0;...
    255,255,51];
ClusterColors = ClusterColors/255;

for iROI = 1:numel(ROI)
    
    clear X_roi
    
    X_roi = X(ROI(iROI).VertOfInt{hs},:);
    
    idx = [];
    C = {};
    sumd = [];
    
    mapping = zeros(NbVertex,1);
    
    for NbCluster=1:MaxNbCluster
        
        fprintf('  Running clustering with %i clusters\n', NbCluster)
        
        [idx(:,NbCluster),C{NbCluster},sumd] = kmeans(X_roi,NbCluster, ...
            'EmptyAction', 'drop', 'Replicates', Replicates, 'Options',opts);

        mapping(ROI(iROI).VertOfInt{hs},NbCluster) = idx(:,NbCluster);
        
        TotalDist(NbCluster) = sum(sumd);
        
    end
    

    %%
    figure('Name', 'Scree plot - Profiles', 'Position', FigDim, ...
        'Color', [1 1 1], 'Visible', Visibility);
    
    subplot(121)
    
    plot(1:MaxNbCluster, TotalDist(1:end), 'linewidth', 2)
    
    NbCluster = find(abs(diff(sign(diff(TotalDist,2))))>1)+2;
    NbCluster = NbCluster(1);
    
    hold on
    plot(NbCluster,TotalDist(NbCluster),'ro');
    
    t=ylabel('Total distance');
    set(t,'fontsize', FontSize);
    
    t=xlabel('Number of clusters');
    set(t,'fontsize', FontSize);
    
    set(gca,'tickdir', 'out', ...
        'xtick', 1:2:length(TotalDist), ...
        'xticklabel', 1:2:MaxNbCluster, ...
        'ticklength', [0.01 0.01], 'fontsize', FontSize)
    
    C = C{NbCluster};

    subplot(122)
    for iCluster=1:NbCluster
        hold on
        plot(C(iCluster,:),'color', 'k')
%         plot(C(iCluster,:),'color', ClusterColors(iCluster,:))
    end
    
    [~,NAME,~] =  fileparts([Name{1} ROI(iROI).name Name{2}]);
    
    mtit(strrep(NAME,'_', ' '),'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, [Name{1} ROI(iROI).name Name{2} '.tif'], '-dtiff')

    write_vtk([Name{1} ROI(iROI).name Name{2} '_cluster_' num2str(NbCluster) '.vtk'], Vertex, Face, mapping(:,NbCluster))

    clear t H TotalDist
    
end
end