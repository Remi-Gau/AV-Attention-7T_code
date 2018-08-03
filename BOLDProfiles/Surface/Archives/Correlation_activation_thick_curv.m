%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

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
%     '15';...
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
DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers-2,1) DesMat'];
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
    
%     GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
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
        
        thickness_file=spm_select('FPList', fullfile(SubjectFolder, 'Structural','CBS'), ...
            ['T1_' SubjID '_.*_' HsSufix 'cr_gm_avg_thickness.vtk']);
        [VertexThick,FaceThick,Thickness{hs}] = read_vtk(thickness_file, 0, 1);
        
        curve_file=spm_select('FPList', fullfile(SubjectFolder, 'Structural','CBS'), ...
            ['T1_' SubjID '_.*_' HsSufix 'cr_gm_avg_curv.vtk']);
        [VertexCurv,FaceCurv,Curvature{hs}] = read_vtk(curve_file, 0, 1);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData; %#ok<*SAGROW>
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
        
    end
    
    
    if any(NbVertex ~= NbVertices)
        %         NbVertex
        %         NbVertices
        error('The number of vertices does not match.')
    end
    
    close all
    
    
    %% Basic condition
    Cdt =[1 1;2 2;3 3];
    ToPlot={'Cst','Lin'};
    
    for iToPlot = 1:numel(ToPlot)
        
        for iCdt = 1:3
            
            X_lh = Curvature{1};
            X_rh = Curvature{2};
            
            Z_lh = Thickness{1};
            Z_rh = Thickness{2};
            
            Y_lh = nan(1,NbVertex(1));
            Y_lh(1,VertexWithDataHS{1}) = BetaCdtion{1}(iToPlot,:,Cdt(iCdt,2));
            Y_rh = nan(1,NbVertex(2));
            Y_rh(1,VertexWithDataHS{2}) = BetaCdtion{2}(iToPlot,:,Cdt(iCdt,2));
            
            for iROI = 1:4%numel(ROI)
                
                X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
                Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
                Z = [Z_lh(ROI(iROI).VertOfInt{1}) Z_rh(ROI(iROI).VertOfInt{2})];

                X(isnan(Y)) =[];
                Z(isnan(Y)) =[];
                Y(isnan(Y)) =[];
                
                R=corrcoef(X,Y);
                rho_curv(iROI,iCdt,SubjInd,iToPlot) = R(1,2);
                beta = glmfit(X, Y, 'normal');
                slope_curv(iROI,iCdt,SubjInd,iToPlot) = beta(2);
                
                R=corrcoef(Z,Y);
                rho_thick(iROI,iCdt,SubjInd,iToPlot) = R(1,2);
                beta = glmfit(Z, Y, 'normal');
                slope_thick(iROI,iCdt,SubjInd,iToPlot) = beta(2);
                
                clear X Y Z

            end
            
            clear Y_lh Y_rh X_lh Y_rh Z_lh Z_rh
            
        end

        
    end
    
    clear Cdt Name ToPlot Range Cdt Name ToPlot  BetaCdtion
    
    
end

if str2double(MatlabVer(1:4))>2013
    if KillGcpOnExit
        delete(gcp)
    end
else
end

save(fullfile(StartFolder,'Results','Profiles','Surfaces','CorrelationActDeactCurvThick.mat'), ...
    'slope_thick', 'rho_thick', 'rho_curv', 'slope_curv')

cd(StartFolder)


return

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

load(fullfile(StartFolder,'Results','Profiles','Surfaces','CorrelationActDeactCurvThick.mat'))

Subj2Include = true(13,1);
Subj2Include([4 11 12]) = false;

Xpos =  [1 1.5];


%%
close all
CdtionName={'thick vs A','thick vs V'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:4
        
        subplot(4,2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
         tmp = squeeze(slope_thick(iROI,iCdt,:,:));

%         tmp = atanh(tmp);
        [H,P] = ttest(tmp);
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),.75,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end

        h = errorbar(Xpos, mean(tmp,1), nanstd(tmp,1), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp, ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -2 1])
        
%         t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(iROI).name));
        t=ylabel(sprintf('%s \n slope',ROI(iROI).name));
        set(t,'fontsize',8)
        
    end
end

subplot(4,2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)

subplot(4,2,2)
t=title(CdtionName{2});
set(t,'fontsize',12)




CdtionName={'thick vs A','thick vs V'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:4
        
        subplot(4,2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
         tmp = squeeze(rho_thick(iROI,iCdt,:,:));

        tmp = atanh(tmp);
        [H,P] = ttest(tmp);
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),.4,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end

        h = errorbar(Xpos, mean(tmp,1), nanstd(tmp,1), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp, ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -.5 .5])
        
        t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(iROI).name));
%         t=ylabel(sprintf('%s \n slope',ROI(iROI).name));
        set(t,'fontsize',8)
        
    end
end

subplot(4,2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)

subplot(4,2,2)
t=title(CdtionName{2});
set(t,'fontsize',12)
 


%%
CdtionName={'curv vs A','curv vs V'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:4
        
        subplot(4,2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
         tmp = squeeze(slope_curv(iROI,iCdt,:,:));

%         tmp = atanh(tmp);
        [H,P] = ttest(tmp);
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),.75,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end

        h = errorbar(Xpos, mean(tmp,1), nanstd(tmp,1), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp, ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -2 1])
        
%         t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(iROI).name));
        t=ylabel(sprintf('%s \n slope',ROI(iROI).name));
        set(t,'fontsize',8)
        
    end
end

subplot(4,2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)

subplot(4,2,2)
t=title(CdtionName{2});
set(t,'fontsize',12)




CdtionName={'curv vs A','curv vs V'};
figure('name', 'Basic', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility)
for iCdt=1:numel(CdtionName)
    for iROI=1:1:4
        
        subplot(4,2,iCdt+2*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
         tmp = squeeze(rho_curv(iROI,iCdt,:,:));

        tmp = atanh(tmp);
        [H,P] = ttest(tmp);
        for i=1:numel(P)
            if P(i)<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P(i));
            end
            t = text(Xpos(i),.4,sprintf(Sig));
            set(t,'fontsize',8);
            
            if H(i)==1
                set(t,'color','r');
            end
        end

        h = errorbar(Xpos, mean(tmp,1), nanstd(tmp,1), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp, ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -.5 .5])
        
        t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(iROI).name));
%         t=ylabel(sprintf('%s \n slope',ROI(iROI).name));
        set(t,'fontsize',8)
        
    end
end

subplot(4,2,1)
t=title(CdtionName{1});
set(t,'fontsize',12)

subplot(4,2,2)
t=title(CdtionName{2});
set(t,'fontsize',12)