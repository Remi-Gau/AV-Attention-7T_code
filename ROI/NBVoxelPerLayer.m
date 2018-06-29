NbLayers = 6;

StartFolder=fullfile(pwd, '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

Include = ones(size(SubjectList,1),1);
% Include([4 11]) = 0;
Include = logical(Include);


% Color for Subjects
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

% COLOR_Subject = repmat([128 128 128],13,1);
COLOR_Subject=COLOR_Subject/255;


FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles');
cd(FigureFolder)
load(fullfile(FigureFolder, 'Volumes', [num2str(NbLayers) '_layers'], ...
    strcat('Data_Median_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')

%%
SavedTxt = 'SizeROI.csv';
fid = fopen (fullfile(StartFolder,SavedTxt), 'w');
for iROI=1:length(AllSubjects_Data)
    fprintf(fid, ',%s,,,', AllSubjects_Data(iROI).name);
end
fprintf(fid, '\n');
for iROI=1:length(AllSubjects_Data)
    fprintf(fid, ',Mean,Std,(Min-,Max)');
end
fprintf(fid, '\n');
for iLayer=2:(NbLayers+1)
    fprintf(fid, 'Layer %i,', iLayer-1);
    for iROI=1:length(AllSubjects_Data)
        fprintf(fid, '%i,', mean(AllSubjects_Data(iROI).VoxelCount(Include,iLayer)));
        fprintf(fid, '%i,', std(AllSubjects_Data(iROI).VoxelCount(Include,iLayer)));
        fprintf(fid, '%i,', min(AllSubjects_Data(iROI).VoxelCount(Include,iLayer)));
        fprintf(fid, '%i,', max(AllSubjects_Data(iROI).VoxelCount(Include,iLayer)));
    end
    fprintf(fid, '\n');
end
fprintf(fid, 'Percent ROI coverage,');
for iROI=1:length(AllSubjects_Data)
    tmp = 100*AllSubjects_Data(iROI).ROI_Coverage(Include,1)./AllSubjects_Data(iROI).ROI_Coverage(Include,2);
    fprintf(fid, '%i,', mean(tmp));
    fprintf(fid, '%i,,,', std(tmp));
end
fclose (fid);


%%
Scatter = linspace(0,.3,11);

Fontsize = 12;

figure('Name', 'Test', 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', 'on');

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for iROI=1:length(AllSubjects_Data)
    
    subplot(2,3,iROI)
 
    hold on
        
    DATA = fliplr(AllSubjects_Data(iROI).VoxelCount(Include,2:end));
       
    for iLayer = 1:size(DATA,2)
        
        h = errorbar(iLayer,mean(DATA(:,iLayer)),std(DATA(:,iLayer)), '.k');
        set(h, 'linewidth', 1, 'MarkerSize', 15)
        
        for SubjInd=1:size(DATA,1)
            plot(iLayer+.2+Scatter(SubjInd), DATA(SubjInd,iLayer), ...
                'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
                'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
        end
    end
    
    title(strrep(AllSubjects_Data(iROI).name, '_', '-'), 'fontsize', Fontsize)
    
    set(gca,'tickdir', 'out', 'xtick', 1:size(DATA,2), 'xticklabel', linspace(0,100,size(DATA,2)), ...
        'ticklength', [0.05 0.05], 'fontsize', Fontsize)
    t=xlabel('Cortical depth');
    set(t,'fontsize',Fontsize);
    t=ylabel('Number of voxels');
    set(t,'fontsize',Fontsize);
    
    h = axis;
    axis([.8 size(DATA,2)+.8 h(3) h(4)])
    
end

print(gcf, fullfile(StartFolder, 'Figures', ['NbVoxels_' num2str(size(DATA,2)) 'Layers.tif']), '-dtiff')
