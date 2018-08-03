function Profile_Analysis_Surf_06_Basic_Plot
% Compare the different ROIs

clear; clc; close all;

Visibility='on';

NbLayers = 6;

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

StartFolder = fullfile(pwd, '..', '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ResultsFolder = fullfile(StartFolder,'Figures','6_layers');


 
%%
clear ROI
ROI(1) = struct('name', 'A1');
ROI(end+1) = struct('name', 'A1_V_act');
ROI(end+1) = struct('name', 'A1_V_deact');
ROI(end+1) = struct('name', 'A1_A_act');
ROI(end+1) = struct('name', 'A1_A_deact');
ROI(end+1) = struct('name', 'A1_act');
ROI(end+1) = struct('name', 'A1_deact');

ROI(end+1) = struct('name', 'PT');
ROI(end+1) = struct('name', 'PT_V_act');
ROI(end+1) = struct('name', 'PT_V_deact');
ROI(end+1) = struct('name', 'PT_A_act');
ROI(end+1) = struct('name', 'PT_A_deact');
ROI(end+1) = struct('name', 'PT_act');
ROI(end+1) = struct('name', 'PT_deact');

ROI(end+1) = struct('name', 'V1');
ROI(end+1) = struct('name', 'V1_V_act');
ROI(end+1) = struct('name', 'V1_V_deact');
ROI(end+1) = struct('name', 'V1_A_act');
ROI(end+1) = struct('name', 'V1_A_deact');
ROI(end+1) = struct('name', 'V1_act');
ROI(end+1) = struct('name', 'V1_deact');

ROI(end+1) = struct('name', 'V2-3');
ROI(end+1) = struct('name', 'V2-3_V_act');
ROI(end+1) = struct('name', 'V2-3_V_deact');
ROI(end+1) = struct('name', 'V2-3_A_act');
ROI(end+1) = struct('name', 'V2-3_A_deact');
ROI(end+1) = struct('name', 'V23_act');
ROI(end+1) = struct('name', 'V23_deact');

All_BOLD_Profile_surf = nan(NbLayers+2,6,numel(ROI),size(SubjectList,1));
All_BOLD_Profile_surf_STD = nan(NbLayers+2,6,numel(ROI),size(SubjectList,1));

for SubjInd = 1:size(SubjectList,1)
    
    
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Surfaces');
    
%     load(fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID],'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')

    
    for ROI_Ind = 1:numel(ROI)
        
        File2Load = fullfile(SubjectFolder,...
            strcat('Data_Surf_Block_', ROI(ROI_Ind).name, '_', num2str(NbLayers+2), '_layers.mat'));
        
        if exist(File2Load, 'file')
            load(File2Load, 'Data_ROI')
            
            All_BOLD_Profile_surf(:,:,ROI_Ind,SubjInd) =  Data_ROI.MEAN;
            All_BOLD_Profile_surf_STD(:,:,ROI_Ind,SubjInd) = squeeze(nanmean(Data_ROI.LayerSEM,2));
        else
            error('File %s is missing.', File2Load)
            All_BOLD_Profile_surf(:,:,ROI_Ind,SubjInd) = nan(NbLayers,6);
            All_BOLD_Profile_surf_STD(:,:,ROI_Ind,SubjInd) = nan(NbLayers,6);
        end
        
    end
end

%%
SubPlots = {[1 2 6 7], 3, 8, 4, 9, 5, 10};

for iROI=1:7:numel(ROI)

    
    %%
    figure('name', [ROI(iROI).name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
        'Color', [1 1 1], 'Visible', Visibility);
    
    tmp2 = squeeze(All_BOLD_Profile_surf(:,:,iROI:iROI+4,:));%-squeeze(All_BOLD_Profile_surf_STD(:,:,iROI:iROI+4,:));
    MIN = min(tmp2(:));
    tmp2 = squeeze(All_BOLD_Profile_surf(:,:,iROI:iROI+4,:));%+squeeze(All_BOLD_Profile_surf_STD(:,:,iROI:iROI+4,:));
    MAX = max(tmp2(:));

    for i=iROI:iROI+6
        
        LEGEND={ROI(i).name}';
        
        subplot(2,5,SubPlots{i-(iROI)+1})
        
        ROIs_DATA.DATA_surf = squeeze(All_BOLD_Profile_surf(:,:,i,:));
        ROIs_DATA.DATA_surf_STD = squeeze(All_BOLD_Profile_surf_STD(:,:,i,:));
        ROIs_DATA.name = ROI(i).name;
        ROIs_DATA.MAX =MAX;
        ROIs_DATA.MIN =MIN;
        
        PlotROIs(ROIs_DATA, LEGEND)
        
        if i==iROI
            legend({...
                'A stim - A att';...
                'V stim - A att';...
                'AV stim - A att';...
                'A stim - V att';...
                'V stim - V att';...
                'AV stim - V att';...
                })
        end
    end

    print(gcf, fullfile(ResultsFolder, [ROI(iROI).name '_surf_all_cdt.tif']), '-dtiff')
    
    
    %%
    Transparent = 'on';
    Fontsize = 10;
    Color = 'rgbrgb';
    
    for i=iROI:iROI+6
        
        ROIs_DATA.DATA_surf = squeeze(All_BOLD_Profile_surf(:,:,i,:));
        MEAN_surf=nanmean(ROIs_DATA.DATA_surf,3);
        SEM_surf=nanstd(ROIs_DATA.DATA_surf,3);
        
        tmp2 = ROIs_DATA.DATA_surf;
        MIN = min(tmp2(:));
        tmp2 = ROIs_DATA.DATA_surf;
        MAX = max(tmp2(:));
        
        figure('name', [ROI(i).name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
            'Color', [1 1 1], 'Visible', Visibility);
        
        iSubplot = 1;
        
        for icdt=1:size(MEAN_surf,2)
            
            subplot(2,3,iSubplot)
            hold on
            grid on
            
            if icdt > 3
                LineStyle = '--';
            else
                LineStyle = '-';
            end
            
            shadedErrorBar(0:size(MEAN_surf,1)-1, flipud(MEAN_surf(:,icdt)), flipud(SEM_surf(:,icdt)), ...
                {'LineStyle', LineStyle, 'LineWidth', 2, 'Color', Color(icdt)}, Transparent)
            
            for iSubj=1:size(ROIs_DATA.DATA_surf,3)
                plot(0:size(MEAN_surf,1)-1, flipud(ROIs_DATA.DATA_surf(:,icdt,iSubj)), 'color',[0.5 0.5 0.5])
            end
            
            plot([0 size(ROIs_DATA.DATA_surf,1)],[0 0], ':k','LineWidth', 1)
            
            axis([-0.1 size(ROIs_DATA.DATA_surf,1)-.9 MIN MAX])
            
            t=ylabel('Beta values');
            set(t,'fontsize',Fontsize); clear t
            
            t=xlabel('% Cortical depth');
            set(t,'fontsize',Fontsize); clear t
            
            set(gca,'tickdir', 'out', 'xtick', 0:size(ROIs_DATA.DATA_surf,1)-1 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_surf,1)))), ...
                'ticklength', [0.01 0.01], 'fontsize', Fontsize-2)
            
            iSubplot=iSubplot+1;
        end
        
        mtit(strrep(ROI(i).name, '_', '-'))
        
        print(gcf, fullfile(ResultsFolder, [ROI(i).name '_surf.tif']), '-dtiff')
    end
end


end



function PlotROIs(ROIs_DATA, LEGEND)

Fontsize = 8;

MEAN_surf=nanmean(ROIs_DATA.DATA_surf,3);
% SEM_surf=nanstd(ROIs_DATA.DATA_surf,3);

MIN = ROIs_DATA.MIN;
MAX = ROIs_DATA.MAX;

Color = 'rgbrgb';


%%
for icdt=1:size(MEAN_surf,2)

    hold on
    grid on
    
    if icdt > 3 
        LineStyle = '--';
    else
        LineStyle = '-';
    end
    
    plot(0:size(MEAN_surf,1)-1,flipud(MEAN_surf(:,icdt)), 'LineStyle', LineStyle, 'LineWidth', 1, 'Color', Color(icdt))
%     shadedErrorBar(0:size(MEAN_surf,1)-1, flipud(MEAN_surf(:,icdt)), flipud(SEM_surf(:,icdt)), ...
%         {'LineStyle', LineStyle, 'LineWidth', 1, 'Color', Color(icdt)}, Transparent)

end

    plot([0 size(ROIs_DATA.DATA_surf,1)],[0 0], ':k','LineWidth', 1)
    
    axis([-0.1 size(ROIs_DATA.DATA_surf,1)-.9 MIN MAX])
        
    t=ylabel('Beta values');
    set(t,'fontsize',Fontsize); clear t
    
    t=xlabel('% Cortical depth');
    set(t,'fontsize',Fontsize); clear t
    
    set(gca,'tickdir', 'out', 'xtick', 0:size(ROIs_DATA.DATA_surf,1)-1 ,'xticklabel', fliplr(round(linspace(0,100,size(ROIs_DATA.DATA_surf,1)))), ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize-2)
    
    title(strrep(LEGEND{1}, '_', '-'))

end