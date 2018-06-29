function Struct_qT1_surf_Plot
% Compare the different ROIs

clear; clc; close all;

Visibility='on';

NbLayers = 11;

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

ROI(1) = struct('name', 'A1');
ROI(end+1) = struct('name', 'PT');
ROI(end+1) = struct('name', 'V1');
ROI(end+1) = struct('name', 'V2-3');
ROI(end+1) = struct('name', 'Cortex');

All_qT1_Profile = nan(NbLayers+1,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    for ROI_Ind = 1:numel(ROI)
        File2Load = fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(ROI_Ind).name, '_surf_', num2str(NbLayers), '.mat'));
        if exist(File2Load, 'file')
            load(File2Load, 'qT1_Profile')
            All_qT1_Profile(:,ROI_Ind,SubjInd) = qT1_Profile(:,1);
        else
            warning('missing data for %s', File2Load)
            All_qT1_Profile(:,ROI_Ind,SubjInd) = nan(NbLayers,1);
        end
        clear qT1_Profile
    end
end


%%
Name = 'Comparison A1 - PT';

LEGEND={ROI(1:2).name}';

ROIs_DATA.DATA = All_qT1_Profile(:,[1:2 end],:);
ROIs_DATA.name = {ROI(1:2).name}';

PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, Visibility)


%%
Name = 'Comparison V1 - V2-3';

LEGEND={ROI(3:4).name}';

ROIs_DATA.DATA = All_qT1_Profile(:,[3:4 end],:);
ROIs_DATA.name = {ROI(3:4).name}';

PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, Visibility)

end


function PlotComparison(ROIs_DATA, Name, LEGEND, SubjectList, Visibility)

Fontsize = 8;

COLOR = [
    0 0 1; ...
    1 0 0; ...
    0 1 0; ...
    0 0 0;];

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

mn = size(ROIs_DATA.DATA,3);
n  = round(mn^0.4);
m  = ceil(mn/n);

MEAN=nanmean(ROIs_DATA.DATA,3);
STD=nanstd(ROIs_DATA.DATA,3);
SEM=nanstd(ROIs_DATA.DATA,3);

Transparent=1;

h = figure('name', [Name '- SUBJECTS'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

g = figure('name', [Name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
    'Color', [1 1 1], 'Visible', Visibility);

for ROI_Ind=1:size(LEGEND,1)
    
    temp = find(strcmp(LEGEND{ROI_Ind}, ROIs_DATA.name));
    
    figure(g)
    
    subplot(1,size(LEGEND,1)+1,ROI_Ind)
    hold on
    shadedErrorBar(0:size(MEAN,1)-1, MEAN(:,temp), SEM(:,temp), ...
        {'LineWidth', 3, 'Color', COLOR(ROI_Ind,:)}, Transparent)
    
    
    subplot(1,size(LEGEND,1)+1,size(LEGEND,1)+1)
    hold on
    grid on
    %     errorbar([1:size(MEAN,1)]+.1*(ROI_Ind-1), MEAN(:,temp), SEM(:,temp), ...
    %         'color', COLOR(ROI_Ind,:), 'LineWidth', 2);
    
    shadedErrorBar(0:size(MEAN,1)-1, MEAN(:,temp), SEM(:,temp), ...
        {'LineWidth', 2, 'Color', COLOR(ROI_Ind,:)}, Transparent);

    
    for SubjInd=1:mn
        
        % Plot a given ROI on the same subplot
        figure(g)
        set(g, 'Visible', Visibility)
        
        subplot(1,size(LEGEND,1)+1,ROI_Ind)
        
        t=title(strrep(LEGEND{ROI_Ind},'_',' '));
        set(t,'fontsize', Fontsize); clear t
        
        hold on
        grid on
        
        %         plot(ROIs_DATA.DATA(:,temp,SubjInd), 'color', COLOR_Subject(SubjInd,:), 'LineWidth', 1)
        plot(0:size(MEAN,1)-1,ROIs_DATA.DATA(:,temp,SubjInd), 'color', [.5 .5 .5], 'LineWidth', 1)
        
        
        
        
        % Plot a given Subject on the same subplot
        figure(h)
        set(h, 'Visible', Visibility)
        subplot(m,n,SubjInd)
        
        t=title(['Subject ' SubjectList(SubjInd,:)]);
        set(t,'fontsize', Fontsize); clear t
        
        hold on
        grid on
        
        plot(0:size(MEAN,1)-1, ROIs_DATA.DATA(:,temp,SubjInd), 'color', COLOR(ROI_Ind,:), 'LineWidth', 2)
        
        axis([-0.5 size(ROIs_DATA.DATA,1)-.5 min(ROIs_DATA.DATA(:)) max(ROIs_DATA.DATA(:))])
        
        t=ylabel('T1 values');
        set(t,'fontsize',Fontsize); clear t
        
        t=xlabel('% Cortical depth');
        set(t,'fontsize',Fontsize); clear t
        
    set(gca,'tickdir', 'out', 'xtick', 1:size(ROIs_DATA.DATA,1)-2 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA,1)-2))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)
        
    end

    clear temp
end


figure(g)

for ROI_Ind=1:size(LEGEND,1)

    subplot(1,size(LEGEND,1)+1,ROI_Ind)
    axis([-0.5 size(ROIs_DATA.DATA,1)-.5 min(ROIs_DATA.DATA(:)) max(ROIs_DATA.DATA(:))])
    
    t=ylabel('T1 values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 1:size(ROIs_DATA.DATA,1)-2 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA,1)-2))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
end

subplot(1,size(LEGEND,1)+1,size(LEGEND,1)+1)

    shadedErrorBar(0:size(MEAN,1)-1, MEAN(:,end), SEM(:,end), ...
        {'LineWidth', 2, 'Color', 'k'}, Transparent);

axis([-0.5 size(ROIs_DATA.DATA,1)-.5 min(ROIs_DATA.DATA(:)) max(ROIs_DATA.DATA(:))])

t=ylabel('T1 values');
set(t,'fontsize',Fontsize); clear t

t=xlabel('% Cortical depth');
set(t,'fontsize',Fontsize); clear t

    set(gca,'tickdir', 'out', 'xtick', 1:size(ROIs_DATA.DATA,1)-2 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA,1)-2))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize)


end