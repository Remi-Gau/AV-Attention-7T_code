%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '14';...
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

DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
% DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);

Bins = 250;
FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

%%
for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
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
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
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
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        Features = mean(FeaturesAll,4);
        
        X=repmat(DesMat,size(Features,3),1);
        
        Y = shiftdim(Features,1);
        Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
        
        BetaAtt{hs}(:,:,end+1) = pinv(X)*Y;
        %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
        
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
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            clear Features Beta2Sel B X Y Mapping
            
        end
        
        clear CondInd AttCondInd CrossSensCondInd CrossSensCondNames AttCondNames BlockInd CondInd
        
    end
    
    
    if any(NbVertex ~= NbVertices)
        %         NbVertex
        %         NbVertices
        error('The number of vertices does not match.')
    end
    
    close all
    
    
    %% Basic condition
    Cdt =[1 2;1 3;2 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1%:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1%:3
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1%:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
%                 subplot(numel(ROI),3,iCdt+3*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpBetaCdt(:,iCdt,iROI,iToPlot, SubjInd) =  Beta;
                
                t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
                set(t,'fontsize',14)
                t=ylabel([Name{Cdt(iCdt,2)} ' stim']);
                set(t,'fontsize',14)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        set(gca,'fontsize', 14)
        
        mtit(['Subject ' SubjID ' - Baseline - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
%         print(gcf, fullfile(FigureFolder,'Baseline', ...
%             ['Subj_' SubjID '_CorrBaseline_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    
    %% Attention condition
    Cdt =[1 2;1 3;2 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1:3
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
                subplot(numel(ROI),3,iCdt+3*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpBetaAtt(:, iCdt,iROI,iToPlot,SubjInd) =  Beta;
                
                t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
                set(t,'fontsize',10)
                t=ylabel([Name{Cdt(iCdt,2)} ' stim']);
                set(t,'fontsize',10)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        mtit(['Subject ' SubjID ' - Attention (A-V) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['Subj_' SubjID '_CorrAtt_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    
    %% Cross sensory to condition
    Cdt =[1 2;2 1];
    Name = {...
        'AV-V','A'; ...
        'AV-A','V'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1:2
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
                subplot(numel(ROI),2,iCdt+2*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpBetaCrossSens(:,iCdt,iROI,iToPlot, SubjInd) =  Beta;
                
                t=xlabel([Name{iCdt,2} ' stim']);
                set(t,'fontsize',10)
                t=ylabel([Name{iCdt,1}]);
                set(t,'fontsize',10)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        mtit(['Subject ' SubjID ' - CrossSensory VS Stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['Subj_' SubjID '_CorrCrossSensory_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    
    %% Attention to stim
    Cdt =[1 4;2 4;3 4];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1:3
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
                subplot(numel(ROI),3,iCdt+3*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpBetaStimAtt(:, iCdt,iROI,iToPlot,SubjInd) =  Beta;
                
                t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
                set(t,'fontsize',10)
                t=ylabel('Att_A - Att_V');
                set(t,'fontsize',10)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        mtit(['Subject ' SubjID ' - Stim vs Attention (A-V) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['Subj_' SubjID '_CorrStimVSAtt_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    
    %% Stim specific attention to stim
    Cdt =[1 1;2 2;3 3];
    Name={'A','V','AV'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1:3
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
                subplot(numel(ROI),3,iCdt+3*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpBetaStimSpeAtt(:, iCdt,iROI,iToPlot,SubjInd) =  Beta;
                
                t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
                set(t,'fontsize',10)
                t=ylabel(sprintf('Att_A-Att_V\nfor %s stim', Name{Cdt(iCdt,1)}));
                set(t,'fontsize',10)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        mtit(['Subject ' SubjID ' - Stim vs Attention - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['Subj_' SubjID '_CorrStimSpeVSAtt_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    
    %% CrossMod to attention
    Cdt =[1 4;2 4];
    Name = {'AV-A'; 'AV-V'};
    ToPlot={'Cst','Lin'};
    Range = [-10 10; -5 5];
    
    for iToPlot = 1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 1:2
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaAtt{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaAtt{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
                subplot(numel(ROI),2,iCdt+2*(iROI-1))
                
                hold on
                Beta = PlotScatterDensity(X,Y,Range(iToPlot,:),Range(iToPlot,:),Bins);
                GrpCrossModAtt(:, iCdt,iROI,iToPlot,SubjInd) =  Beta;
                
                t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
                set(t,'fontsize',10)
                t=ylabel('Att_A - Att_V');
                set(t,'fontsize',10)
                title(ROI(iROI).name)
                
                clear X Y
                
            end
            
            clear Y_lh Y_rh X_lh Y_rh
            
        end
        
        mtit(['Subject ' SubjID ' - CrossMod vs Attention - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['Subj_' SubjID '_CorrCrossModVSAtt_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range Cdt Name ToPlot Range BetaCdtion BetaCrossSens BetaAtt
    
    
end

if str2double(MatlabVer(1:4))>2013
    if KillGcpOnExit
        delete(gcp)
    end
else
end

save(fullfile(StartFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'), ...
    'ROI', 'GrpBetaStimSpeAtt', 'GrpBetaStimAtt', 'GrpCrossModAtt', 'GrpBetaCrossSens','GrpBetaAtt','GrpBetaCdt')

cd(StartFolder)


%%
StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

load(fullfile(StartFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'))

Subj2Include = true(13,1);
Subj2Include([4 11]) = false;

Xpos =  [1 1.5];


%%
close all
CdtionName={'A vs V','A vs AV','V vs AV'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCdt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),1.75,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.5 1.75])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',12)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder,'Baseline', 'Grp_CorrBaseline.tif'), '-dtiff')


%%
close all
CdtionName={'A vs AV-V','V vs AV-A'};
figure('name', 'CrossSensory', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCrossSens(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),1.5,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.25 1.5])
    end
end

subplot(numel(ROI),2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),2,3)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,5)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,7)
t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder,'CrossSensory', 'Grp__CorrCrossSensory.tif'), '-dtiff')


%%
close all
CdtionName={'Att Mod_{A stim} vs Att Mod_{V stim}','Att Mod_{A stim} vs Att Mod_{AV stim}','Att Mod_{V stim} vs AttMod_{AV stim}'};
figure('name', 'Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.55,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.20 0.55])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_CorrAttention.tif'), '-dtiff')


%% Stim vs Att
close all
CdtionName={'A stim vs Att Mod','V stim vs Att Mod','AV Stim vs AttMod'};
figure('name', 'Stiim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaStimAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_Corr_Stim_VS_Attention.tif'), '-dtiff')


%% Stim vs stim specific Att
close all
CdtionName={'A stim vs Att Mod A stim','V stim vs Att Mod V stim','AV Stim vs AttMod AV stim'};
figure('name', 'Stim VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),3,iCdt+3*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaStimSpeAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),3,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),3,3)
t=title(CdtionName{3});
set(t,'fontsize',8)

subplot(numel(ROI),3,4)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),3,10)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'Attention', 'Grp_Corr_StimSpe_VS_Attention.tif'), '-dtiff')



%% Cross vs Att
    Name = {'AV-A'; 'AV-V'};
close all
CdtionName={'AV-A vs Att Mod','AV-V vs Att Mod'};
figure('name', 'CrossMod VS Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI),2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpCrossModAtt(3,iCdt,iROI,:,Subj2Include));
        tmp = atanh(tmp);
        [H,P] = ttest(tmp');
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),0.25,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end
        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -1:.25:1 ,'yticklabel', -1:.25:1,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.75 0.25])
    end
end

subplot(numel(ROI),2,1)
t=title(CdtionName{1});
set(t,'fontsize',8)
t=ylabel(sprintf('%s \n Corr Coef (Fishcer trans)',ROI(1).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,2)
t=title(CdtionName{2});
set(t,'fontsize',8)

subplot(numel(ROI),2,3)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(2).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,5)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(3).name));
set(t,'fontsize',8)

subplot(numel(ROI),2,7)
t=ylabel(sprintf('%s \nCorr Coef\n(Fishcer trans)',ROI(4).name));
set(t,'fontsize',8)

print(gcf, fullfile(FigureFolder, 'CrossSensory', 'Grp_Corr_CrossMod_VS_Attention.tif'), '-dtiff')
