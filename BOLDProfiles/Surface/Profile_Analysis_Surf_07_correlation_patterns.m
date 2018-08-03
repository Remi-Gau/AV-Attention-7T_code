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

DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
% DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);


for SubjInd = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder,'BetaMapping','8Surf');
    
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
    fprintf(' Averaging for ROI:\n')
    
    for iROI = 1:2 %numel(ROI)
        
        clear Data_ROI
        
        Data_ROI.name = ROI(iROI).name;
        
        fprintf(['  '  Data_ROI.name '\n'])
        
        Features = cat(1,Features_lh(ROI(iROI).VertOfInt{1},:,:), ...
            Features_rh(ROI(iROI).VertOfInt{2},:,:));
        
        FeaturesL = Features_lh(ROI(iROI).VertOfInt{1},:,:);
        FeaturesR = Features_rh(ROI(iROI).VertOfInt{2},:,:);
        
        % (A, AV) x 2 (Att V, Att A)
        % 1. A AttA and A Att V
        % 2. AV AttA and AV Att V
        % should be bigger than for
        % 1.  A AttA and AV AttA
        % 2. A Att V  and AV Att V
        Correlation = [...
            1,4;
            3,6;
            1,3;
            4,6];
        
        for iCor = 1:size(Correlation,1)
            
            Beta2Sel = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{Correlation(iCor,1)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            
            tmp = Features(:,2:7,Beta2Sel);
            tmp(any(any(isnan(tmp),2),3),:,:) = [];
            
            X=repmat(DesMat,size(tmp,3),1);
            Y = shiftdim(tmp,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            B1 = pinv(X)*Y;
            

            Beta2Sel = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{Correlation(iCor,2)} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            
            tmp2 = Features(:,2:7,Beta2Sel);
            tmp2(any(any(isnan(tmp2),2),3),:,:) = [];
            
            X=repmat(DesMat,size(tmp2,3),1);
            Y = shiftdim(tmp2,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            B2 = pinv(X)*Y;
            
            
            RHO_cst(SubjInd,iCor,iROI) = corr(B1(1,:)',B2(1,:)','type','Spearman');
            RHO_lin(SubjInd,iCor,iROI) = corr(B1(2,:)',B2(2,:)','type','Spearman');
            
            
            tmp = mean(tmp,3);
            tmp2 = mean(tmp2,3);

            for iLayer=1:6
                RHO_Layer(SubjInd,iCor,iLayer,iROI) = corr(tmp(:,iLayer),tmp2(:,iLayer),'type','Spearman');
            end

            RHO(SubjInd,iCor,iROI) = corr(tmp(:),tmp2(:),'type','Spearman');
            
        end

    end
    
    cd(Results_Folder)
    %         save(strcat('Data_Surf_Block_', Data_ROI.name, '_', num2str(NbLayers), '_layers.mat'), 'Data_ROI')
    
end

return

%%
close all

Dim = [0.5 4.5 0.15 .7];
FigDim = [50 50 1600 650];
XAxis = {'A_A / A_V'
     'AV_A / AV_V'
     'A_A / AV_A'
     'A_V / AV_V'};

SubPlot=1;

figure('name', 'RHO', 'position', FigDim)

for iROI=1:2
    
subplot(2,4,SubPlot)
errorbar(1:4,mean(RHO(:,:,iROI)),nansem(RHO(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:4, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(1:4,mean(RHO_cst(:,:,iROI)),nansem(RHO_cst(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:4, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(1:4,mean(RHO_lin(:,:,iROI)),nansem(RHO_lin(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:4, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(repmat(1:6,4,1)'+repmat([0:.1:.3]',1,6)',...
    squeeze(mean(RHO_Layer(:,:,:,iROI)))',...
    squeeze(nansem(RHO_Layer(:,:,:,iROI)))', 'marker', '.' , 'markersize', 10)
set(gca, 'xtick', 1:6, 'xticklabel', 1:6, 'fontsize', 8)
legend(XAxis, 'Location', 'northwest')
axis([0.5 6.5 0.15 .65])
SubPlot=SubPlot+1;


end

subplot(2,4,1)
title('Whole ROI')
ylabel('A1')
subplot(2,4,5)
ylabel('PT')
subplot(2,4,2)
title('Cst')
subplot(2,4,3)
title('Lin')
subplot(2,4,4)
title('Layers')



%%
close all

RHO_f = [mean(atanh(RHO(:,1:2,:)),2) mean(atanh(RHO(:,3:4,:)),2)];
RHO_cst_f = [mean(atanh(RHO_cst(:,1:2,:)),2) mean(atanh(RHO_cst(:,3:4,:)),2)];
RHO_lin_f = [mean(atanh(RHO_lin(:,1:2,:)),2) mean(atanh(RHO_lin(:,3:4,:)),2)];
RHO_Layer_f = [mean(atanh(RHO_Layer(:,1:2,:,:)),2) mean(atanh(RHO_Layer(:,3:4,:,:)),2)];

Dim = [0.5 2.5 0.3 .65];
FigDim = [50 50 1600 650];
XAxis = {'(A+AV)_A vs (A+AV)_V'
     'A_{A+V} vs AV_{A+V}'};

SubPlot=1;

figure('name', 'RHO', 'position', FigDim)

for iROI=1:2
    
subplot(2,4,SubPlot)
errorbar(1:2,mean(RHO_f(:,:,iROI)),nansem(RHO_f(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:2, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(1:2,mean(RHO_cst_f(:,:,iROI)),nansem(RHO_cst_f(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:2, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(1:2,mean(RHO_lin_f(:,:,iROI)),nansem(RHO_lin_f(:,:,iROI)), '.k' , 'markersize', 10)
set(gca, 'xtick', 1:2, 'xticklabel', XAxis, 'fontsize', 8)
axis(Dim)
SubPlot=SubPlot+1;

subplot(2,4,SubPlot)
errorbar(repmat(1:6,2,1)'+repmat([0 .1]',1,6)',...
    squeeze(mean(RHO_Layer_f(:,:,:,iROI)))',...
    squeeze(nansem(RHO_Layer_f(:,:,:,iROI)))', 'marker', '.' , 'markersize', 10)
set(gca, 'xtick', 1:6, 'xticklabel', 1:6, 'fontsize', 8)
legend(XAxis, 'Location', 'northwest')
axis([0.5 6.5 0.15 .75])
SubPlot=SubPlot+1;


end

subplot(2,4,1)
title('Whole ROI')
ylabel('A1')
subplot(2,4,5)
ylabel('PT')
subplot(2,4,2)
title('Cst')
subplot(2,4,3)
title('Lin')
subplot(2,4,4)
title('Layers')

%% run permutation test on A1 cst
sets = {};
for iSub=1:11
    sets{iSub} = [-1 1]; %#ok<*AGROW>
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


tmp = diff(RHO_cst_f(:,:,1),[],2);
for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:);
    tmp2 = repmat(tmp2',1,size(tmp,2));
    Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
end
P = sum( ...
    abs( Perms ) > ...
    repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
    / size(Perms,1)

