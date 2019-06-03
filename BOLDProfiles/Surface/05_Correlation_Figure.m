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


DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
% DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);

Bins = 250;
FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

%%
for SubjInd = 2 %1:size(SubjectList,1)
    
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
    
    


    %% Cross sensory to condition
    close all
    
    Cdt =[1 2;2 1];
    Name = {...
        'AV-V','A'; ...
        'AV-A','V'};
    ToPlot={'Cst','Lin'};
    Range = [-3 3; -3 3];
    
    for iToPlot = 2 %1:numel(ToPlot)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
        
        for iCdt = 2
            
            X_lh = nan(1,NbVertex(1));
            X_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,1));
            X_rh = nan(1,NbVertex(2));
            X_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,1));
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaCrossSens{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaCrossSens{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1 %:numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                
%                 subplot(numel(ROI),2,iCdt+2*(iROI-1))
                
                hold on
                
%                 scatter(X,Y,'.')
                Beta = PlotScatterDensity2(X,Y,Range(iToPlot,:),Range(iToPlot,:),100);
%                 GrpBetaCrossSens(:,iCdt,iROI,iToPlot, SubjInd) =  Beta;
                
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
        
%         print(gcf, fullfile(FigureFolder,'CrossSensory', ...
%             ['Subj_' SubjID '_CorrCrossSensory_' ToPlot{iToPlot} '.tif']), '-dtiff')
        
    end
    
    clear Cdt Name ToPlot Range
    
    


    
%     clear Cdt Name ToPlot Range Cdt Name ToPlot Range BetaCdtion BetaCrossSens BetaAtt
    
    
end

