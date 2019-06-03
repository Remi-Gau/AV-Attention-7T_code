%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

Print=0;

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
Ind = NbLayers:-1:1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

FigDim = [100 100 1800 1000];
Visible = 'off';

NbWorkers = 1;

% MatlabVer = version('-release');
% if str2double(MatlabVer(1:4))>2013
%     if isempty(gcp)
%         KillGcpOnExit = 1;
%         parpool(NbWorkers);
%     else
%         KillGcpOnExit = 0;
%     end
%     % else
%     %     if matlabpool('size') == 0
%     %         KillGcpOnExit = 1;
%     %         matlabpool(NbWorkers)
%     %     elseif matlabpool('size') ~= NbWorkers
%     %         matlabpool close
%     %         matlabpool(NbWorkers)
%     %         KillGcpOnExit = 0;
%     %     else
%     %         KillGcpOnExit = 0;
%     %     end
% end


for SubjInd = 3 %1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
%     GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
%         ['Subject_' SubjID], 'FFX_Block');
    
Data_Folder = fullfile(SubjectFolder,'BetaMapping','8Surf');
%     Data_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
%         ['Subject_' SubjID],'BetaMapping','8Surf');
    
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
        
        InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
            ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
        
%         [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);
%         NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData;
            MappingBothHS{hs} = AllMapping;
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
            
            
            save(FeatureSaveFile,'Vertex','Face','AllMapping','VertexWithData', '-v7.3')
            
            VertexWithDataHS{hs} = VertexWithData;
            MappingBothHS{hs} = AllMapping;
        end
        
        clear Betas InfSurfFile Mapping A B C
        
    end
    
%     if any(NbVertex ~= NbVertices)
%         NbVertex
%         NbVertices
%         error('The number of vertices does not match.')
%     end
    
    
    Features_lh = nan(NbVertex(1),NbLayers,size(MappingBothHS{1},3));
    % Features_lh = [];
    Features_lh(VertexWithDataHS{1},:,:) = MappingBothHS{1};
    
    Features_rh = nan(NbVertex(2),NbLayers,size(MappingBothHS{2},3));
    % Features_rh = [];
    Features_rh(VertexWithDataHS{2},:,:) = MappingBothHS{2};
    
    %%
    fprintf(' Averaging for ROI:\n')
    
    for iROI = 1:numel(ROI)
        
        clear Data_ROI
        
        Data_ROI.name = ROI(iROI).name;
        
        fprintf(['  '  Data_ROI.name '\n'])
        
        Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
            Features_rh(ROI(iROI).VertOfInt{2},:,:));
        
        FeaturesL = Features_lh(ROI(iROI).VertOfInt{1},:,:);
        FeaturesR = Features_rh(ROI(iROI).VertOfInt{2},:,:);
        
        if Print
            figure('position', FigDim, 'name', 'Attention', 'Color', [1 1 1], 'visible', Visible)
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            
            set(gcf, 'PaperUnits','centimeters');
            set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        end
        
        for CondInd = 4%:length(Conditions_Names) % For each Condition
            
            Beta2Sel = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            
            if Print
                % Plot distribution
                for LayerInd = 1:NbLayers
                    DistToPlot{LayerInd} = nanmean(Features(:,LayerInd,Beta2Sel),3); %#ok<SAGROW>
                end
                
                subplot(2,3,CondInd)
                hold on
                distributionPlot(DistToPlot, 'xValues', 1:numel(DistToPlot), 'color', 'k', ...
                    'distWidth', 0.8, 'showMM', 1, 'globalNorm', 2, 'histOpt', 1.1)
                set(gca,'tickdir', 'out', 'xtick', 1:numel(DistToPlot) , ...
                    'xticklabel', Ind, ...
                    'ticklength', [0.01 0.01], 'fontsize', 12)
                plot([0 numel(DistToPlot)+1], [0 0], '--k')
                axis([0 numel(DistToPlot)+1 -10 10])
            end
            
            for BlockInd = 1:numel(Beta2Sel) % For each Block
                
                tmp = Features(:,:,Beta2Sel(BlockInd));
                tmpL = FeaturesL(:,:,Beta2Sel(BlockInd));
                tmpR = FeaturesR(:,:,Beta2Sel(BlockInd));
                
                if CondInd==1 && BlockInd==1
                    fprintf('  NaNs: %i ; Zeros: %i\n',sum(any(isnan(tmp),2)),sum(any(tmp==0,2)))
                    Data_ROI.NaNorZero = [any(isnan(tmp),2) any(tmp==0,2)];
                end
                
                %                 tmp(any([any(isnan(tmp),2) any(tmp==0,2)],2),:) = [];
                %                 tmpL(any([any(isnan(tmpL),2) any(tmpL==0,2)],2),:) = [];
                %                 tmpR(any([any(isnan(tmpR),2) any(tmpR==0,2)],2),:) = [];
                
                tmp(any(tmp==0,2),:) = [];
                tmpL(any(tmpL==0,2),:) = [];
                tmpR(any(tmpR==0,2),:) = [];
                
                for LayerInd = 1:NbLayers % Averages over voxels of a given layer
                    
                    Data_ROI.LayerMean(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmean(tmp(:,LayerInd));
                    Data_ROI.LayerMedian(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmedian(tmp(:,LayerInd));
                    Data_ROI.LayerSTD(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanstd(tmp(:,LayerInd));
                    Data_ROI.LayerSEM(Ind(LayerInd),BlockInd,CondInd)  = ...
                        nansem(tmp(:,LayerInd));
                    
                    Data_ROI.LayerMeanL(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmean(tmpL(:,LayerInd));
                    Data_ROI.LayerMedianL(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmedian(tmpL(:,LayerInd));
                    Data_ROI.LayerSTD_L(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanstd(tmpL(:,LayerInd));
                    Data_ROI.LayerSEM_L(Ind(LayerInd),BlockInd,CondInd)  = ...
                        nansem(tmpL(:,LayerInd));
                    
                    Data_ROI.LayerMeanR(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmean(tmpR(:,LayerInd));
                    Data_ROI.LayerMedianR(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanmedian(tmpR(:,LayerInd));
                    Data_ROI.LayerSTD_R(Ind(LayerInd),BlockInd,CondInd) = ...
                        nanstd(tmpR(:,LayerInd));
                    Data_ROI.LayerSEM_R(Ind(LayerInd),BlockInd,CondInd)  = ...
                        nansem(tmpR(:,LayerInd));
                    
                end
                
            end
            
        end
        
        Data_ROI.MEAN=squeeze(nanmean(Data_ROI.LayerMean,2));
        Data_ROI.MEDIAN=squeeze(nanmean(Data_ROI.LayerMedian,2));
        Data_ROI.STD=squeeze(nanstd(Data_ROI.LayerMean,2));
        Data_ROI.SEM=squeeze(nansem(Data_ROI.LayerMean,2));
        
        Data_ROI.MEAN_L=squeeze(nanmean(Data_ROI.LayerMeanL,2));
        Data_ROI.MEDIAN_L=squeeze(nanmean(Data_ROI.LayerMedianL,2));
        Data_ROI.STD_L=squeeze(nanstd(Data_ROI.LayerMeanL,2));
        Data_ROI.SEM_L=squeeze(nansem(Data_ROI.LayerMeanL,2));
        
        Data_ROI.MEAN_R=squeeze(nanmean(Data_ROI.LayerMeanR,2));
        Data_ROI.MEDIAN_R=squeeze(nanmean(Data_ROI.LayerMedianR,2));
        Data_ROI.STD_R=squeeze(nanstd(Data_ROI.LayerMeanR,2));
        Data_ROI.SEM_R=squeeze(nansem(Data_ROI.LayerMeanR,2));
        
        if Print
            % figure
            subplot(2,3,1)
            t=ylabel('A attention');
            set(t,'fontsize',12);
            t=title('A stimulation');
            set(t,'fontsize',12);
            
            subplot(2,3,2)
            t=title('V stimulation');
            set(t,'fontsize',12);
            
            subplot(2,3,3)
            t=title('AV stimulation');
            set(t,'fontsize',12);
            
            subplot(2,3,4)
            t=ylabel('V attention');
            set(t,'fontsize',12);
            
            
            %         mtit(['Surf-Subj ' SubjID '-' strrep(Data_ROI.name,'_','-')],'xoff',0,'yoff',.025)
            print(gcf, fullfile(Results_Folder, ...
                ['Suibj_' SubjID '_VoxelsDist_surf_', Data_ROI.name, '_', num2str(NbLayers), '_layers.tif']), '-dtiff')
        end
        
        cd(Results_Folder)
        save(strcat('Data_Surf_Block_', Data_ROI.name, '_', num2str(NbLayers), '_layers.mat'), 'Data_ROI')
        
    end
    
end

if str2double(MatlabVer(1:4))>2013
    if KillGcpOnExit
        delete(gcp)
    end
else
end

