function Struct_qT1_vol_surf_comp
% Compare the different ROIs

clear; clc; close all;

Visibility='on';

NbLayers = 10;

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

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder, 'Figures', 'Structural', 'T1_profiles');


Suffix = ' - Both HS ';

%%
ROI(1) = struct('name', 'A1');
ROI(end+1) = struct('name', 'PT_thres');

ROI(end+1) = struct('name', 'V1_thres');
ROI(end+1) = struct('name', 'V2-3_thres');

ROI(end+1) = struct('name', 'Cortex');

All_qT1_Profile_vol = nan(NbLayers,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    for ROI_Ind = 1:numel(ROI)
        File2Load = fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(ROI_Ind).name, '_vol_', num2str(NbLayers), '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'qT1_Profile')
            All_qT1_Profile_vol(:,ROI_Ind,SubjInd) = qT1_Profile(:,1);
        else
            error('File %s is missing.', File2Load)
            All_qT1_Profile_vol(:,ROI_Ind,SubjInd) = nan(NbLayers,1);
        end
        
    end
end

All_qT1_Profile_GLM = nan(NbLayers,numel(ROI)-1,size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    for ROI_Ind = 1:numel(ROI)-1
        File2Load = fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(ROI_Ind).name, '_GLM_', num2str(NbLayers), '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'qT1_Profile')
            All_qT1_Profile_GLM(:,ROI_Ind,SubjInd) = qT1_Profile(:,1);
        else
            error('File %s is missing.', File2Load)
            All_qT1_Profile_GLM(:,ROI_Ind,SubjInd) = nan(NbLayers,1);
        end
        
    end
end


%%
clear ROI
ROI(1) = struct('name', 'A1');
ROI(end+1) = struct('name', 'PT');

ROI(end+1) = struct('name', 'V1');
ROI(end+1) = struct('name', 'V2-3');

ROI(end+1) = struct('name', 'Cortex');

All_qT1_Profile_surf = nan(NbLayers+2,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    for ROI_Ind = 1:numel(ROI)
        File2Load = fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(ROI_Ind).name, '_surf_', num2str(NbLayers+1), '.mat'));
        %             strcat('qT1_profile_', ROI(ROI_Ind).name, '_rh_surf_', num2str(NbLayers+1), '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'qT1_Profile')
            All_qT1_Profile_surf(:,ROI_Ind,SubjInd) = qT1_Profile(:,1);
        else
            error('File %s is missing.', File2Load)
            All_qT1_Profile_surf(:,ROI_Ind,SubjInd) = nan(NbLayers,1);
        end
        
    end
end





%%
Name = ['Comparison A1 - PT' Suffix];

LEGEND={ROI([1:2 end]).name}';

ROIs_DATA.DATA_vol = All_qT1_Profile_vol(:,[1:2 end],:);
ROIs_DATA.DATA_surf = All_qT1_Profile_surf(:,[1:2 end],:);
ROIs_DATA.DATA_GLM = All_qT1_Profile_GLM(:,[1:2 end],:);
ROIs_DATA.name = {ROI([1:2 end]).name}';

PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, FigureFolder, Visibility)


%%
Name = ['Comparison V1 - V2-3' Suffix];

LEGEND={ROI([3:4 end]).name}';

ROIs_DATA.DATA_vol = All_qT1_Profile_vol(:,[3:4 end],:);
ROIs_DATA.DATA_surf = All_qT1_Profile_surf(:,[3:4 end],:);
ROIs_DATA.DATA_GLM = All_qT1_Profile_GLM(:,[3:4 end],:);
ROIs_DATA.name = {ROI([3:4 end]).name}';

PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, FigureFolder, Visibility)

return

end



function PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, FigureFolder, Visibility)

Fontsize = 8;

COLOR = [
    0 0 1; ...
    1 0 0; ...
    0 0 0];

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
COLOR_Subject=COLOR_Subject/255;

mn = size(ROIs_DATA.DATA_vol,3);
n  = round(mn^0.4);
m  = ceil(mn/n);

MEAN_vol=nanmean(ROIs_DATA.DATA_vol,3);
STD_vol=nanstd(ROIs_DATA.DATA_vol,3);
SEM_vol=nanstd(ROIs_DATA.DATA_vol,3);

MEAN_GLM=nanmean(ROIs_DATA.DATA_GLM,3);
STD_GLM=nanstd(ROIs_DATA.DATA_GLM,3);
SEM_GLM=nanstd(ROIs_DATA.DATA_GLM,3);

MEAN_surf=nanmean(ROIs_DATA.DATA_surf,3);
STD_surf=nanstd(ROIs_DATA.DATA_surf,3);
SEM_surf=nanstd(ROIs_DATA.DATA_surf,3);

Transparent=1;


%%
g = figure('name', [Name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for ROI_Ind=1:size(LEGEND,1)
    
    temp = find(strcmp(LEGEND{ROI_Ind}, ROIs_DATA.name));
    
    subplot(1,size(LEGEND,1),ROI_Ind)
    hold on
    shadedErrorBar(1:size(MEAN_vol,1), MEAN_vol(:,temp), SEM_vol(:,temp), ...
        {'LineWidth', 1, 'Color', COLOR(ROI_Ind,:)}, Transparent)
    if ROI_Ind<3
        shadedErrorBar(1:size(MEAN_GLM,1), MEAN_GLM(:,temp), SEM_GLM(:,temp), ...
            {'LineStyle', ':', 'LineWidth', 1, 'Color', COLOR(ROI_Ind,:)}, Transparent)
    end
    shadedErrorBar(0:size(MEAN_surf,1)-1, MEAN_surf(:,temp), SEM_surf(:,temp), ...
        {'LineStyle', '--', 'LineWidth', 1, 'Color', COLOR(ROI_Ind,:)}, Transparent)
    
    
    t=title(strrep(LEGEND{ROI_Ind},'_',' '));
    set(t,'fontsize', Fontsize); clear t
    
    
    axis([0.5 size(ROIs_DATA.DATA_vol,1)+.5 min(ROIs_DATA.DATA_vol(:)) max(ROIs_DATA.DATA_vol(:))])
    
    t=ylabel('T1 values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 1:size(ROIs_DATA.DATA_vol,1) ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_vol,1)))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    axis([0.5 size(ROIs_DATA.DATA_vol,1)+.5 1400 2600])
    
end

mtit(Name, 'xoff',0,'yoff',.025)

print(gcf, fullfile(FigureFolder, strrep([Name '.tif'],' ','')), '-dtiff')


%%
h = figure('name', [Name '- SUBJECTS'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for SubjInd=1:mn
    
    
    % Plot a given Subject on the same subplot
    
    subplot(m,n,SubjInd)
    
    t=title(['Subject ' SubjectList(SubjInd,:)]);
    set(t,'fontsize', Fontsize); clear t
    
    hold on
    grid on
    
    for ROI_Ind=1:size(LEGEND,1)-1
        temp = find(strcmp(LEGEND{ROI_Ind}, ROIs_DATA.name));
        plot(ROIs_DATA.DATA_vol(:,temp,SubjInd), 'color', COLOR(ROI_Ind,:), 'LineWidth', 1)
        plot(ROIs_DATA.DATA_GLM(:,temp,SubjInd), 'color', COLOR(ROI_Ind,:), 'LineWidth', 1,  'LineStyle', ':')
        plot(0:size(MEAN_surf,1)-1,ROIs_DATA.DATA_surf(:,temp,SubjInd), 'color', COLOR(ROI_Ind,:), 'LineWidth', 1,  'LineStyle', '--')
    end
    
    t=ylabel('T1 values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 1:size(ROIs_DATA.DATA_vol,1) ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_vol,1)))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    axis([0.5 size(ROIs_DATA.DATA_vol,1)+.5 1400 2600])
    
end

mtit([Name '- Subjects'], 'xoff',0,'yoff',.025)

print(gcf, fullfile(FigureFolder, strrep([Name '_SUBJECTS.tif'],' ','')), '-dtiff')

end