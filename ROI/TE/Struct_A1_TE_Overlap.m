clear; clc; close all;

StartFolder=pwd;
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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
    %     '15';...
    '16'
    ];

MinLayer = 1;
NbLayers = 6;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

Mask.ROI(1) = struct('name', 'A1', 'fname', 'A1.nii');
Mask.ROI(end+1) = struct('name', 'TE1.0', 'fname', 'rrwTE_1.0_MNI.nii');
Mask.ROI(end+1) = struct('name', 'TE1.1', 'fname', 'rrwTE_1.1_MNI.nii');
% Mask.ROI(end+1) = struct('name', 'TE1.2', 'fname', 'rrwTE_1.2_MNI.nii');


%% Get overlap between ROIs
for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n  Processing subject : ' SubjID '\n'])
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    Data_Folder = fullfile(SubjectFolder, 'Transfer');
    ROI_Folder = fullfile(Data_Folder, 'ROI');
    
    %%
    LayerVol = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', ...
        [TargetLayerFile '.nii']);
    LayerVolHdr = spm_vol(LayerVol);
    Cortex = logical(spm_read_vols(spm_vol(LayerVolHdr)));
    
    %%
    cd(ROI_Folder)
    ROI_Vols = spm_read_vols(spm_vol(char({Mask.ROI.fname}')));
    
    for iImg=1:size(ROI_Vols,4)
        tmp = Cortex;
        tmp(:,:,:,2) = ROI_Vols(:,:,:,iImg);
        ROI_Vols(:,:,:,iImg) = all(tmp,4);
    end
    
    %%
    for iImg=1:size(ROI_Vols,4)
        Results(SubjInd,iImg) = numel(find(ROI_Vols(:,:,:,iImg))); %#ok<*SAGROW>
    end
    
    Results(SubjInd,4) = numel(find(all(ROI_Vols(:,:,:,1:2),4)));
    Results(SubjInd,5) = numel(find(all(ROI_Vols(:,:,:,[1 3]),4)));
    Results(SubjInd,6) = numel(find(all(ROI_Vols(:,:,:,[2 3]),4)));
    Results(SubjInd,7) = numel(find(all(ROI_Vols(:,:,:,[1:3]),4)));
    
    Results %#ok<NOPTS>
    
end

cd(StartFolder)

% save('VoxelOverlap.mat', 'Results')


%% Plot
COLOR_Subject= [
    0,0,0;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

Fontsize = 8;

Visible = 'on';

mn = size(SubjectList,1);
n  = round(mn^0.4);
m  = ceil(mn/n);

load('VoxelOverlap.mat', 'Results')

cd(fullfile(StartFolder, 'Figures'))


%% Venn Diagram
figure('name', 'Overlap', 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visible);

COLOR = {'r','g','b'};
   

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf([' Plotting subject : ' SubjID '\n'])

    subplot(n,m,SubjInd)
    t=title(['Subject ' SubjectList(SubjInd,:)]);
    set(t,'fontsize', Fontsize); clear t
    hold on
    [~, S] = venn(Results(SubjInd,:), 'ErrMinMode','ChowRodgers', 'FaceColor', COLOR);
    
    %Now label each zone
    for i = 1:length(Results(SubjInd,:))
        t=text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(sprintf('%2.1f',Results(SubjInd,i)/10000)));
        set(t,'fontsize', Fontsize); clear t
    end
    clear i
    
    set(gca,'tickdir', 'out', 'xtick', [], 'xticklabel', [], ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    axis off
    axis image
    
end



print(gcf, 'Overlap_A1_TE1.0_TE1.1.tif', '-dtiffnocompression')


%%
TMP = nanmean(Results);

figure('name', 'Overlap', 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visible);
[~, S] = venn(TMP, 'ErrMinMode','ChowRodgers', 'FaceColor', {'r','g','b'});
axis off

for i = 1:7
    t=text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), sprintf('%2.1f', TMP(i)/10000));
    set(t,'fontsize', Fontsize); clear t
end

axis off
axis image

LEGEND = char({Mask.ROI.name}');

legend('Location','NorthEastOutside', LEGEND)
legend boxoff

print(gcf, 'MeanOverlap_A1_TE1.0_TE1.1.tif', '-dtiffnocompression')


