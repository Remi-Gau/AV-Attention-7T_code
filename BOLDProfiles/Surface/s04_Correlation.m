% Computes vertex wise correlation for different contrasts.
% Gives a scatter plot of those correlations for each subject
% Correlation are not cross-validated across days

%%
clc; clear;

CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))

FigureFolder = fullfile(CodeFolder,'Figures');

SourceFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

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

ROIs_to_use = 1:4;


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

opt.Bins = 250;
opt.FigDim = [100, 100, 1000, 1500];
opt.Visibility = 'off';
opt.Range = [-10 10; -5 5];
opt.FigureFolder = FigureFolder;
opt.fontsize = 14;
opt.print = 0;


%%
GrpBetaStimSpeAtt = [];
GrpBetaStimAtt = [];
GrpCrossModAtt = [];
GrpBetaCrossSens = [];
GrpBetaAtt = [];
GrpBetaCdt = [];

for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(SourceFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
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
    load(fullfile(SubjectFolder,'ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    ROI = ROI(ROIs_to_use);
    
    %% Read features
    fprintf(' Reading VTKs\n')
    
    % Format for reading the vertices from the VTK file
    Spec = repmat('%f ', 1, NbLayers);
    
    NbVertices = nan(1,2);
    
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
        [AllMapping, Face, Vertex, VertexWithData] = ...
            load_vtk_beta_maps(FeatureSaveFile, Data_Folder, HsSufix, NbVertices, NbLayers);
        
        VertexWithDataHS{hs} = VertexWithData;
        
        clear Vertex Face
        
        %% Run GLMs for A, V and AV
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
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaCdtion{hs}(:,:,CondInd) = pinv(X)*Y;
            
            clear Features Beta2Sel Beta2Sel2 X Y Mapping
            
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
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaAtt{hs}(:,:,AttCondInd) = pinv(X)*Y;
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        Features = mean(FeaturesAll,4);
        
        X=repmat(DesMat,size(Features,3),1);
        
        Y = shiftdim(Features,1);
        Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
        
        BetaAtt{hs}(:,:,end+1) = pinv(X)*Y;
        
        clear FeaturesAll
        
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
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            BetaCrossSens{hs}(:,:,CrossSensCondInd) = pinv(X)*Y;
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        clear CondInd AttCondInd CrossSensCondInd CrossSensCondNames AttCondNames BlockInd CondInd
        
    end
    
    close all
    
    %% Basic condition
    opt.Cdt = [1 2; 1 3; 2 3];
    opt.Name={'A','V', 'AV'};
    opt.ToPlot={'Cst', 'Lin'};
    
    GrpBetaCdt = plot_vertex_wise_correlation(BetaCdtion, BetaCdtion, GrpBetaCdt, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    %% Attention condition
    opt.Cdt = [1 2; 1 3; 2 3];
    opt.Name={'A_{att_A-att_V}','V_{att_A-att_V}','AV_{att_A-att_V}'};
    
    GrpBetaAtt = plot_vertex_wise_correlation(BetaAtt, BetaAtt, GrpBetaAtt, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    
    %% Cross sensory to condition
    opt.Cdt =[1 2; 2 1];
    opt.Name = {...
        'AV-V','A'; ...
        'AV-A','V'};
    
    GrpBetaCrossSens = plot_vertex_wise_correlation(BetaCdtion, BetaCrossSens, GrpBetaCrossSens, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    
    %% Overall effect of attention to stim
    opt.Cdt =[1 4; 2 4 ;3 4];
    opt.Name={'A','V','AV', 'all'};
    
    GrpBetaStimAtt = plot_vertex_wise_correlation(BetaCdtion, BetaAtt, GrpBetaStimAtt, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    
    %% Stim specific attention to stim
    opt.Cdt =[1 1; 2 2; 3 3];
    opt.Name={'A','V','AV'};
    
    GrpBetaStimSpeAtt = plot_vertex_wise_correlation(BetaCdtion, BetaAtt, GrpBetaStimSpeAtt, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    
    %% CrossMod to attention
    opt.Cdt =[1 4; 2 4];
    opt.Name = {'AV-A'; 'AV-V'; ''; 'All'};
    
    GrpCrossModAtt = plot_vertex_wise_correlation(BetaCrossSens, BetaAtt, GrpCrossModAtt, ...
        NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt);
    
    
    clear BetaCdtion BetaCrossSens BetaAtt VertexWithDataHS NbVertex
    
end

save(fullfile(SourceFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'), ...
    'ROI', 'GrpBetaStimSpeAtt', 'GrpBetaStimAtt', 'GrpCrossModAtt', ...
    'GrpBetaCrossSens','GrpBetaAtt','GrpBetaCdt')


