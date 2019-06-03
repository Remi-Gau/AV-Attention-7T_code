%%
clc; clear;

COLOR =   [...
    255 150 150; ...
    150 255 150; ...
    150 220 255; ...
    255 75 75; ...
    75 255 75; ...
    75 75 255];
COLOR=COLOR/255;

FigDim = [100 100 1800 1000];
Visible = 'off';

MinLayer = 1;
NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

ANTs = 0;

StartFolder=fullfile(pwd, '..', '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
    %     '03';...
    '04';...
    %     '06';...
    %     '07';...
    %     '08';...
    %     '09';...
    %     '11';...
    %     '12';...
    '13';...
    %     '14';...
    '15';...
    %     '16'
    ];



Mask.ROI(1) = struct('name', 'V1_surf_thres', 'fname', 'V1_surf_thres.nii');
Mask.ROI(end+1) = struct('name', 'V2-3_surf_thres', 'fname', 'V2-3_surf_thres.nii');

% Mask.ROI(end+1) = struct('name', 'V1_A_Deact', 'fname', 'V1_surf_thres_A_Stim<Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V1_V_AV_Deact', 'fname', 'V1_surf_thres_V_AV_Stim<Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V1_V_AV_Act', 'fname', 'V1_surf_thres_V_AV_Stim>Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V2-3_A_Deact', 'fname', 'V2-3_surf_thres_A_Stim<Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V2-3_V_AV_Deact', 'fname', 'V2-3_surf_thres_V_AV_Stim<Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V2-3_V_AV_Act', 'fname', 'V2-3_surf_thres_V_AV_Stim>Baseline_thres_0.1.nii');

% Mask.ROI(end+1) = struct('name', 'V1_A_Act', 'fname', 'V1_surf_thres_A_Stim>Baseline_thres_0.1.nii');
% Mask.ROI(end+1) = struct('name', 'V2-3_A_Act', 'fname', 'V2-3_surf_thres_A_Stim>Baseline_thres_0.1.nii');

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

Fontsize = 10;

SkimThres = 20;


for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    %     GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder, 'Transfer');
    
    %     BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
    %         ['Subject_' SubjID],'UpsampledBetas');
    MaskFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'Transfer');
    BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'ManCoreg','V','exp-0000','exp-0000-A','Transform_Volume');
%     cd(BackUpFolder)
%     gunzip('*.gz')
    
    ROI_Folder = fullfile(Data_Folder, 'ROI');
    
    
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Volumes', TargetLayerFile);
    [~,~,~] = mkdir(Results_Folder);
    
    Destination_Folder = fullfile('/home/rxg243/Dropbox/PhD/Experiments/AV_Integration_7T',...
        'Subjects_Data', ['Subject_' SubjID], 'Transfer', 'Profiles', TargetLayerFile);
    
    
    %% Gets ROIs
    
    % Gets global mask from GLM and ROI masks for the data
    Mask.global.hdr = spm_vol(fullfile(MaskFolder, 'rmask.nii'));
    Mask.global.img = logical(spm_read_vols(Mask.global.hdr));
    
    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            Mask.ROI(iROI).hdr = spm_vol(fullfile(ROI_Folder, Mask.ROI(iROI).fname));
        end
    end
    
    hdr = cat(1, Mask.ROI.hdr);
    sts = spm_check_orientations([Mask.global.hdr; hdr]);
    if sts ~= 1
        error('Images not in same space!');
    end
    clear sts hdr i
    
    % Create mask in XYZ format (both world and voxel coordinates)
    [X, Y, Z] = ind2sub(size(Mask.global.img), find(Mask.global.img));
    Mask.global.XYZ = [X'; Y'; Z']; % XYZ format
    clear X Y Z
    Mask.global.size = size(Mask.global.XYZ, 2);
    Mask.global.XYZmm = Mask.global.hdr.mat(1:3,:) ...
        * [Mask.global.XYZ; ones(1, Mask.global.size)]; % voxel to world transformation
    
    % Combine masks
    fprintf(' Open ROI images\n')
    xY.def = 'mask';
    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            xY.spec = fullfile(ROI_Folder, Mask.ROI(iROI).fname);
            [xY, Mask.ROI(iROI).XYZmm, j] = spm_ROI(xY, Mask.global.XYZmm);
            Mask.ROI(iROI).XYZ = Mask.global.XYZ(:,j);
            Mask.ROI(iROI).size = size(Mask.ROI(iROI).XYZ, 2);
            A = spm_read_vols(Mask.ROI(iROI).hdr);
            A(isnan(A)) = 0;
            Mask.ROI(iROI).size(2) = sum(logical(A(:)));
            clear A
        end
    end
    clear xY j i
    
    
    %% Open layer file and gets indices of voxels in each layer
    
    fprintf(' Open layer label images & identify layer label for each voxel\n')
    
    LayerVol = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', ...
        [TargetLayerFile '.nii']);
    
    LayerVolHdr = spm_vol(LayerVol);
    
    % Number of voxel of each ROI and intersection of ROI
    for iROI=1:length(Mask.ROI)
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            Mask.ROI(iROI).Layers = spm_get_data(LayerVolHdr, Mask.ROI(iROI).XYZ);
            for iLayer = 1:max(Mask.ROI(iROI).Layers)
                Mask.ROI(iROI).LayersVox(1,iLayer) = sum(Mask.ROI(iROI).Layers==iLayer);
            end
        end
    end
    
    clear LayerVol LayerVolHdr
    
    %% Opens beta images
    fprintf(' Identifying the relevant beta images\n')
    
    load(fullfile(GLM_Folder, 'SPM.mat'));
    
    BetaNames = char(SPM.xX.name');
    
    BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
        strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
        2);
    
    BetaOfInterest = find(BetaOfInterest);
    
    BetaNames(:,1:6)=[];
    
    % Scaling factor for PSC
    xBF.dt = SPM.xBF.dt;
    xBF.name = SPM.xBF.name;
    xBF.length = SPM.xBF.length;
    xBF.order = SPM.xBF.order;
    xBF = spm_get_bf(xBF);
    SF = max(xBF.bf(:,1));
    
    RegNumbers = GetRegNb(SPM);
    
    clear SPM BetaImg
    
    for CondInd = 1:length(Conditions_Names)
        
        tmp=BetaNames(BetaOfInterest,1:length(Conditions_Names{CondInd}));
        
        BetaImgInd{CondInd}=BetaOfInterest(strcmp(Conditions_Names{CondInd}, cellstr(tmp)));  %#ok<*SAGROW>
        
        for i=1:length(BetaImgInd{CondInd})
            if ANTs
                BetaImg{CondInd}(i,:)=fullfile(BackUpFolder, ['rbeta_' sprintf('%3.4d', BetaImgInd{CondInd}(i)) '_def.nii']);
            else
                BetaImg{CondInd}(i,:)=fullfile(BackUpFolder, ['r4beta_' sprintf('%3.4d', BetaImgInd{CondInd}(i)) '_clone_transform.nii']); %#ok<*UNRCH>
            end
        end
        
        clear tmp i
    end
    clear BetaNames BetaOfInterest
    
    
    %% Opens session constant images
    %     cd(SubjectFolder)
    %     CstFiles = dir('rrSess*');
    %     CstFiles = [repmat([fullfile(SubjectFolder) filesep], [4,1]) char({CstFiles.name}')];
    %     CstFiles = spm_vol(CstFiles);
    
    %% AVERAGING
    
    fprintf(' Averaging for ROI:\n')
    
    load('/data/AV_Integration_2/Figures/Behavioral/BehavioralResults.mat', 'GroupResults')
    load(fullfile('/data','AV_Integration_2','Figures','Behavioral','ExcludeBlocks.mat'), 'ExcludeBlocks')
    
    %     ToKeep = cat(4,GroupResults(SubjInd).FARateBlock(1:3,1:2,:)>0, ...
    %         ExcludeBlocks(SubjInd).ExcludeBlocks, ...
    %         GroupResults(SubjInd).MISS_Block_TOTAL(1:3,1:2,:)>0) ;
    
    ToKeep = cat(4,GroupResults(SubjInd).FARateBlock(1:3,1:2,:)>0, ...
        ExcludeBlocks(SubjInd).ExcludeBlocks) ;
    
    ToKeep = any(ToKeep,4);
    
    for iROI=1:length(Mask.ROI)
        
        if exist(fullfile(ROI_Folder, Mask.ROI(iROI).fname), 'file')
            
            clear Data_ROI
            
            fprintf(['  ' Mask.ROI(iROI).name '\n'])
            
            Data_ROI.name = Mask.ROI(iROI).name;
            Data_ROI.Voxelcount = [Mask.ROI(iROI).size(1)  Mask.ROI(iROI).LayersVox];
            spe = repmat('%i ', [1,max(Mask.ROI(iROI).Layers)+1]);
            fprintf(['   Voxel per layer ' spe '\n'], Data_ROI.Voxelcount)
            
            Data_ROI.ROI_Coverage = Mask.ROI(iROI).size;
            fprintf('   Coverage %0.2f\n', Mask.ROI(iROI).size(1)/Mask.ROI(iROI).size(2))
            
            if ~any(Data_ROI.ROI_Coverage==0) && sum(Data_ROI.Voxelcount~=0)==(NbLayers+1)
                
                
                figure('position', FigDim, 'name', 'Attention', 'Color', [1 1 1], 'visible', Visible)
                set(gca,'units','centimeters')
                pos = get(gca,'Position');
                ti = get(gca,'TightInset');
                
                set(gcf, 'PaperUnits','centimeters');
                set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                set(gcf, 'PaperPositionMode', 'manual');
                set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
                
                
                for CondInd = 1:length(Conditions_Names) % For each Condition
                    
                    clear SourceVolume SourceVolumeStruc
                    
                    SourceVolumeStruc = spm_vol(BetaImg{CondInd});
                    SourceVolume = spm_get_data(SourceVolumeStruc, Mask.ROI(iROI).XYZ);
                    %                     Cst = spm_get_data(CstFiles, Mask.ROI(iROI).XYZ);
                    
                    % Plot distribution
                    for LayerInd = 1:max(Mask.ROI(iROI).Layers)
                        Vox2Sel = Mask.ROI(iROI).Layers==LayerInd;
                        DistToPlot{LayerInd} = mean(SourceVolume(:,Vox2Sel));
                    end
                    subplot(2,3,CondInd)
                    hold on
                    distributionPlot(DistToPlot, 'xValues', 1:max(Mask.ROI(iROI).Layers), 'color', 'k', ...
                        'distWidth', 0.8, 'showMM', 1, 'globalNorm', 2, 'histOpt', 1.1)
                    set(gca,'tickdir', 'out', 'xtick', 1:max(Mask.ROI(iROI).Layers) , ...
                        'xticklabel', max(Mask.ROI(iROI).Layers):-1:1, ...
                        'ticklength', [0.01 0.01], 'fontsize', 12)
                    plot([0 max(Mask.ROI(iROI).Layers)+1], [0 0], '--k')
                    axis([0 max(Mask.ROI(iROI).Layers)+1 -10 10])
                    
                    
                    
                    for BlockInd = 1:size(SourceVolume,1) % For each Block
                        
                        % find the session this block belongs to
                        [Sess,~,~] = find(RegNumbers==BetaImgInd{CondInd}(BlockInd));
                        
                        for LayerInd = 0:max(Mask.ROI(iROI).Layers) % Averages over voxels of a given layer
                            
                            Vox2Sel = Mask.ROI(iROI).Layers==LayerInd;
                            
                            Data_ROI.LayerMean(LayerInd+1,BlockInd,CondInd) = ...
                                nanmean(SourceVolume(BlockInd,Vox2Sel));
                            Data_ROI.LayerMedian(LayerInd+1,BlockInd,CondInd) = ...
                                nanmedian(SourceVolume(BlockInd,Vox2Sel));
                            Data_ROI.LayerSTD(LayerInd+1,BlockInd,CondInd) = ...
                                nanstd(SourceVolume(BlockInd,Vox2Sel));
                            Data_ROI.LayerSEM(LayerInd+1,BlockInd,CondInd)  = ...
                                nansem(SourceVolume(BlockInd,Vox2Sel));
                            
                            % Skim the X percentile of voxels
                            % if the mean is above zero we remove the top X
                            % percent voxels
                            if nanmean(SourceVolume(BlockInd,Vox2Sel))>0
                                Thres = prctile(SourceVolume(BlockInd,Vox2Sel), 100-SkimThres);
                                BottomXVox = all([Vox2Sel;SourceVolume(BlockInd,:)<Thres]);
                            else % Otherwise we remove the X bottom one
                                Thres = prctile(SourceVolume(BlockInd,Vox2Sel), SkimThres);
                                BottomXVox = all([Vox2Sel;SourceVolume(BlockInd,:)>Thres]);
                            end
                            
                            Data_ROI.LayerMeanSkim(LayerInd+1,BlockInd,CondInd) = ...
                                nanmean(SourceVolume(BlockInd,BottomXVox));
                            Data_ROI.LayerMedianSkim(LayerInd+1,BlockInd,CondInd) = ...
                                nanmedian(SourceVolume(BlockInd,BottomXVox));
                            Data_ROI.LayerSTDSkim(LayerInd+1,BlockInd,CondInd) = ...
                                nanstd(SourceVolume(BlockInd,BottomXVox));
                            Data_ROI.LayerSEMSkim(LayerInd+1,BlockInd,CondInd)  = ...
                                nansem(SourceVolume(BlockInd,BottomXVox));
                            
                            Data_ROI.ThresSkim(LayerInd+1,BlockInd,CondInd) = Thres;
                            
                            % Percent signal change
                            %                             PSC = SourceVolume(BlockInd,Vox2Sel).*SF./Cst(Sess,Vox2Sel).*100;
                            %                             Data_ROI.LayerMeanPSC(LayerInd+1,BlockInd,CondInd) = nanmean(PSC);
                            %                             Data_ROI.LayerMedianPSC(LayerInd+1,BlockInd,CondInd) = nanmedian(PSC);
                            %                             Data_ROI.LayerSTDPSC(LayerInd+1,BlockInd,CondInd) = nanstd(PSC);
                            %                             Data_ROI.LayerSEMPSC(LayerInd+1,BlockInd,CondInd) = nansem(PSC);
                            
                        end
                        clear tmp LayerInd
                        
                        % Now do the same over all the layer of the ROI
                        Vox2Sel = Mask.ROI(iROI).Layers>0;
                        
                        Data_ROI.LayerMean(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanmean(SourceVolume(BlockInd,Vox2Sel));
                        Data_ROI.LayerMedian(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanmedian(SourceVolume(BlockInd,Vox2Sel));
                        Data_ROI.LayerSTD(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanstd(SourceVolume(BlockInd,Vox2Sel));
                        Data_ROI.LayerSEM(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd)  = ...
                            nansem(SourceVolume(BlockInd,Vox2Sel));
                        
                        % Skim the X percentile of voxels
                        % if the mean is above zero we remove the top X
                        % percent voxels
                        if nanmean(SourceVolume(BlockInd,Vox2Sel))>0
                            Thres = prctile(SourceVolume(BlockInd,Vox2Sel), 100-SkimThres);
                            BottomXVox = all([Vox2Sel;SourceVolume(BlockInd,:)<Thres]);
                        else % Otherwise we remove the X bottom one
                            Thres = prctile(SourceVolume(BlockInd,Vox2Sel), SkimThres);
                            BottomXVox = all([Vox2Sel;SourceVolume(BlockInd,:)>Thres]);
                        end
                        
                        Data_ROI.LayerMeanSkim(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanmean(SourceVolume(BlockInd,BottomXVox));
                        Data_ROI.LayerMedianSkim(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanmedian(SourceVolume(BlockInd,BottomXVox));
                        Data_ROI.LayerSTDSkim(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = ...
                            nanstd(SourceVolume(BlockInd,BottomXVox));
                        Data_ROI.LayerSEMSkim(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd)  = ...
                            nansem(SourceVolume(BlockInd,BottomXVox));
                        
                        Data_ROI.ThresSkim(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = Thres;
                        
                        % Percent signal change
                        %                         PSC = SourceVolume(BlockInd,Vox2Sel).*SF./Cst(Sess,Vox2Sel).*100;
                        %                         Data_ROI.LayerMeanPSC(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = nanmean(PSC);
                        %                         Data_ROI.LayerMedianPSC(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = nanmedian(PSC);
                        %                         Data_ROI.LayerSTDPSC(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = nanstd(PSC);
                        %                         Data_ROI.LayerSEMPSC(max(Mask.ROI(iROI).Layers)+2,BlockInd,CondInd) = nansem(PSC);
                        
                    end
                    clear SourceVolume
                    
                    
                    Data_ROI.LayerMeanRestrict(:,:,CondInd) = Data_ROI.LayerMean(:,:,CondInd);
                    Data_ROI.LayerMedianRestrict(:,:,CondInd) = Data_ROI.LayerMedian(:,:,CondInd);
                    [y,x] = ind2sub([3 2],CondInd);
                    Data_ROI.LayerMeanRestrict(:,(squeeze(ToKeep(y,x,:))),CondInd) = NaN;
                    Data_ROI.LayerMedianRestrict(:,(squeeze(ToKeep(y,x,:))),CondInd) = NaN;
                    
                end
                
                
                % Mean over blocks
                Data_ROI.MEAN=squeeze(nanmean(Data_ROI.LayerMean,2));
                Data_ROI.MEDIAN=squeeze(nanmean(Data_ROI.LayerMedian,2));
                Data_ROI.STD=squeeze(nanstd(Data_ROI.LayerMean,2));
                Data_ROI.SEM=squeeze(nansem(Data_ROI.LayerMean,2));
                
                Data_ROI.MEANSkim=squeeze(nanmean(Data_ROI.LayerMeanSkim,2));
                Data_ROI.MEDIANSkim=squeeze(nanmean(Data_ROI.LayerMedianSkim,2));
                Data_ROI.STDSkim=squeeze(nanstd(Data_ROI.LayerMeanSkim,2));
                Data_ROI.SEMSkim=squeeze(nansem(Data_ROI.LayerMeanSkim,2));
                
                %                 Data_ROI.MEANPSC=squeeze(nanmean(Data_ROI.LayerMeanPSC,2));
                %                 Data_ROI.MEDIANPSC=squeeze(nanmean(Data_ROI.LayerMedianPSC,2));
                %                 Data_ROI.STDPSC=squeeze(nanstd(Data_ROI.LayerMeanPSC,2));
                %                 Data_ROI.SEMPSC=squeeze(nansem(Data_ROI.LayerMeanPSC,2));
                
                Data_ROI.MEANRestrict=squeeze(nanmean(Data_ROI.LayerMeanRestrict,2));
                Data_ROI.MEDIANRestrict=squeeze(nanmean(Data_ROI.LayerMedianRestrict,2));
                Data_ROI.STDRestrict=squeeze(nanstd(Data_ROI.LayerMeanRestrict,2));
                Data_ROI.SEMRestrict=squeeze(nansem(Data_ROI.LayerMeanRestrict,2));
                
                
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
                
                if ANTs
                    mtit(['Subj ' SubjID '-' strrep(Data_ROI.name,'_','-'), '- ANTs'],'xoff',0,'yoff',.025)
                    print(gcf, fullfile(Results_Folder, ...
                        ['Subj_' SubjID '_VoxelsDist_', Data_ROI.name, '_', TargetLayerFile, '_ANTs.tif']), '-dtiff')
                else
                    mtit(['Subj ' SubjID '-' strrep(Data_ROI.name,'_','-')],'xoff',0,'yoff',.025)
                    print(gcf, fullfile(Results_Folder, ...
                        ['Subj_' SubjID '_VoxelsDist_', Data_ROI.name, '_', TargetLayerFile '.tif']), '-dtiff')
                end
                
                close all
            end
            
            cd(Results_Folder)
            if ANTs
                save(strcat('Data_Block_', Data_ROI.name, '_', TargetLayerFile, '_ANTs.mat'), 'Data_ROI')
                %                 copyfile(strcat('Data_Block_', Data_ROI.name, '_', TargetLayerFile, '_ANTs.mat'), ...
                %                     Destination_Folder)
            else
                save(strcat('Data_Block_', Data_ROI.name, '_', TargetLayerFile, '.mat'), 'Data_ROI')
                %                 copyfile(strcat('Data_Block_', Data_ROI.name, '_', TargetLayerFile, '.mat'), ...
                %                     Destination_Folder)
            end
            
            clear CondInd
            
        end
        
    end
    
    
    
end


cd(StartFolder)
