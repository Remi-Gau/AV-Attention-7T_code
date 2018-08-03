%%
clc; clear; close all

V1_2_A1 = 1;
Suffix = 'V1_2_A1';

StartDirectory = fullfile(pwd, '..', '..', '..');
cd(StartDirectory)
addpath(genpath(fullfile(StartDirectory, 'SubFun')))
Get_dependencies('D:\Dropbox')

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

FigureFolder = fullfile(StartDirectory, 'Figures', strcat(num2str(NbLayers), '_layers'));

CdtVec = repmat(1:6,12,1);
CdtVec = CdtVec(:)';

for iSub = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(iSub,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder, 'BetaMapping','8Surf');
    
    
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
    load(fullfile(SubjectFolder,'Transfer','ROI',...
        ['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
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
    
    % run a laminar GLM to get the constant of each block and vertices
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
    
    % run a laminar GLM to get the constant for each condition
    Cst_SeedROI = [];
    for iCdt = 1:6
        for i=1:size(tmp,3)
            X=DesMat;
            Y = tmp(:,:,i,iCdt);
            Y=Y';
            B=pinv(X)*Y;
            Cst_SeedROI = [Cst_SeedROI;B(1,:)];
        end
    end
    
    clear tmp Y X i B Features_lh Features_rh Features
    
    %% Plot
    figure('name', 'dist', 'position', [50 50 1200 1200])
    clf
    clc
    X = linspace(-20,20,1000);
    Subplot = 1;
    
    for i=[1 2 4 5]
        
        Row2Select = ismember(CdtVec,i)';
        
        for iROItoPlot=1:2
            
            if iROItoPlot==1
                Data_to_plot = mean(Cst_SeedROI(Row2Select,:));
            else
                Data_to_plot = mean(Cst_TargetROI(Row2Select,:));
            end
            
            subplot(7,2,Subplot)
            hold on
            h = hist(Data_to_plot,X);
            h = bar(X,h);
            set(h, 'FaceColor', 'k', 'EdgeColor', 'k')
            axis tight
            ax = axis;
            
            h = normpdf(X,mean(Data_to_plot),std(Data_to_plot));
            h = h/sum(h)*numel(Data_to_plot);
            h = plot(X,h,'r','linewidth',2);
            
            plot([mean(Data_to_plot) mean(Data_to_plot)], [0 ax(4)], '-b',...
                'linewidth',2)
            plot([median(Data_to_plot) median(Data_to_plot)], [0 ax(4)], '--b',...
                'linewidth',2)
            plot([0 0], [0 ax(4)], 'g')
            
            axis([-10 10 0 ax(4)])
            Subplot=Subplot+1;
        end
        
        
    end
    
    subplot(7,2,1)
    title('V1')
    ylabel('Stim A_{Att A}')
    legend({'Data';'Norm Dis';'Mean';'Median';'Zero'})
    subplot(7,2,2)
    title('A1')
    
    subplot(7,2,3)
    ylabel('Stim V_{Att A}')
    subplot(7,2,5)
    ylabel('Stim A_{Att V}')
    subplot(7,2,7)
    ylabel('Stim V_{Att V}')
    
    Subplot=Subplot+2;
    
    
    for i =1:2
        
        if i==1
            Cdt = [1 2];
        else
            Cdt = [4 5];
        end

        for iROItoPlot=1:2
            
            if iROItoPlot==1
                Data_to_plot = mean(Cst_SeedROI(ismember(CdtVec,Cdt(1))',:))-...
                    mean(Cst_SeedROI(ismember(CdtVec,Cdt(2))',:));
            else
                Data_to_plot = mean(Cst_TargetROI(ismember(CdtVec,Cdt(1))',:))-...
                    mean(Cst_TargetROI(ismember(CdtVec,Cdt(2))',:));
            end
            
            subplot(7,2,Subplot)
            hold on
            h = hist(Data_to_plot,X);
            h = bar(X,h);
            set(h, 'FaceColor', 'k', 'EdgeColor', 'k')
            axis tight
            ax = axis;
            
            h = normpdf(X,mean(Data_to_plot),std(Data_to_plot));
            h = h/sum(h)*numel(Data_to_plot);
            h = plot(X,h,'r','linewidth',2);
            
            plot([mean(Data_to_plot) mean(Data_to_plot)], [0 ax(4)], '-b',...
                'linewidth',2)
            plot([median(Data_to_plot) median(Data_to_plot)], [0 ax(4)], '--b',...
                'linewidth',2)
            plot([0 0], [0 ax(4)], 'g')
            
            axis([-10 10 0 ax(4)])
            Subplot=Subplot+1;
        end
    end
    
    subplot(7,2,11)
    ylabel('(A-V)_{Att A}')
    subplot(7,2,13)
    ylabel('(A-V)_{Att V}')
    
    mtit(['Distriubutions Subject ' SubjID], 'xoff', 0, 'yoff', +0.04 )
    
    print(gcf, fullfile(FigureFolder,['Distriubutions_Subject_' SubjID '.tif']), '-dtiff')
    
    
    
    %% compute connectivity
    fprintf(' Compute connectivity\n')
    % CdtVec: column vector ((72 X 1) to identify which condition on which row
    % Cst_SeedROI: Cst for each vertex of V1 for each block of each condition
    %   (72 X Nb vertices)
    % Cst_TargetROI: Cst for each vertex of A1 for each block of each condition
    %   (72 X Nb vertices)
    
    % Pooled over A and V only separately for each attention condition
    for iCdt = 1:2
        
        % specify the rows to select for the concatenation
        if iCdt==1
            Row2Select = ismember(CdtVec,1:2)'; % A and V stim under A att
        else
            Row2Select = ismember(CdtVec,4:5)'; % A and V stim under V att
        end

        
        % Do regression after averaging
        Y = nanmean(Cst_TargetROI(Row2Select,:),2); % For A1 average across vertices of activation (now a 1 X 24 vector)
        X = nanmean(Cst_SeedROI(Row2Select,:),2); % Same for V1
        B_att_V1toA1(iSub,iCdt,:,1) = glmfit(X-mean(X),Y,'normal'); % Regression V1 to A1
        B_att_A1toV1(iSub,iCdt,:,1) = glmfit(Y-mean(Y),X,'normal'); % Regression A1 to V1

        
        % Do regression before averaging
        % V1 to A1
        for iVert = 1:size(Cst_SeedROI,2) % regress for each vertex
            tmp(:,iVert) = glmfit(...
                Cst_SeedROI(Row2Select,iVert)-nanmean(Cst_SeedROI(Row2Select,iVert)),...
                Y,'normal');
        end
        B_att_V1toA1(iSub,iCdt,:,2)=mean(tmp,2); % average
        clear tmp
        
        % A1 to V1
        for iVert = 1:size(Cst_TargetROI,2) % regress for each vertex
            tmp(:,iVert) = glmfit(...
                Cst_TargetROI(Row2Select,iVert)-nanmean(Cst_TargetROI(Row2Select,iVert)),...
                X,'normal');
        end
        B_att_A1toV1(iSub,iCdt,:,2)=mean(tmp,2); % average
        clear tmp
        
    end

    cd(Results_Folder)
    
    clear Cst_SeedROI Cst_TargetROI iCdt
    
    
end

%%

save(fullfile(Results_Folder,strcat('Connectivity_all_mean_Surf_', num2str(NbLayers), '_layers.mat')), ...
    'B_att_A1toV1', 'B_att_V1toA1')


%%
clc
load(fullfile(pwd,'Connectivity_all_mean_Surf_8_layers.mat'))


%%
fprintf('Connectivity: AVERAGE then REGRESS\n')
fprintf('Connectivity V1 to A1 for A att: %f +/- %f\n', ...
    mean(B_att_V1toA1(:,1,2,1)), std(B_att_V1toA1(:,1,2,1)))
fprintf('Connectivity V1 to A1 for V att: %f +/- %f\n', ...
    mean(B_att_V1toA1(:,2,2,1)), std(B_att_V1toA1(:,2,2,1)))
P = SignPermTest(B_att_V1toA1(:,1,2,1)-B_att_V1toA1(:,2,2,1))

fprintf('\n')
fprintf('Connectivity A1 to V1 for A att: %f +/- %f\n', ...
    mean(B_att_A1toV1(:,1,2,1)), std(B_att_A1toV1(:,1,2,1)))
fprintf('Connectivity A1 to V1 for V att: %f +/- %f\n', ...
    mean(B_att_A1toV1(:,2,2,1)), std(B_att_A1toV1(:,2,2,1)))
P = SignPermTest(B_att_A1toV1(:,1,2,1)-B_att_A1toV1(:,2,2,1))

fprintf('\n\n\n')
fprintf('Connectivity: REGRESS then AVERAGE\n')
fprintf('Connectivity V1 to A1 for A att: %f +/- %f\n', ...
    mean(B_att_V1toA1(:,1,2,2)), std(B_att_V1toA1(:,1,2,2)))
fprintf('Connectivity V1 to A1 for V att: %f +/- %f\n', ...
    mean(B_att_V1toA1(:,2,2,2)), std(B_att_V1toA1(:,2,2,2)))
P = SignPermTest(B_att_V1toA1(:,1,2,2)-B_att_V1toA1(:,2,2,2))

fprintf('\n')
fprintf('Connectivity A1 to V1 for A att: %f +/- %f\n', ...
    mean(B_att_A1toV1(:,1,2,2)), std(B_att_A1toV1(:,1,2,2)))
fprintf('Connectivity A1 to V1 for V att: %f +/- %f\n', ...
    mean(B_att_A1toV1(:,2,2,2)), std(B_att_A1toV1(:,2,2,2)))
P = SignPermTest(B_att_A1toV1(:,1,2,2)-B_att_A1toV1(:,2,2,2))




