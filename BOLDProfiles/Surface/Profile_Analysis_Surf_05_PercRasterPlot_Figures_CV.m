clear
clc

% StartFolder='/home/rxg243/Dropbox/PhD/Experiments/AV_Integration_7T/Archives/';
% StartFolder='/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
StartFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives';
addpath(genpath(fullfile(StartFolder, 'SubFun')))

Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

load(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster_CV.mat'))

Subj2Include = true(11,1);

for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
% ToPermute = [];


% Color map
X = 0:0.001:1;
ColorMap = seismic(size(X,2));


FigDim = [100, 100, 1000, 700];
Visibility = 'on';


Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

Color(10,:) = [1,0,0];

ToPlot={'Cst','Lin'};

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

CorticalDepth = round(linspace(100,0,NbLayers+2));
CorticalDepth([1 end]) = [];


ROI_Groups = {1:4};
ROI_Groups_names = {'A1-PT-V123'};


load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')



%% Grp level  ; raster Att = f(Percentile of V)
close all
mkdir(fullfile(FigureFolder,'Attention'));
NameSwitch={'Att_V-A', 'Astim_Att_V-A', 'Vstim_Att_V-A'};

SubPlot = [3 4 1 2];
iSubPlot = 0;


for iToPlot = 1
    if iToPlot==1
        CLIM = [-.8 .8];
        CLIM2 = CLIM;
    else
        CLIM = [-0.4 .4];
        CLIM2 = CLIM*6;
    end
    
    for iCdt = 3 %1:3
        
        for iROI = 2%:4 % 1:4
            
            figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_A,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_AttfV{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            % plot main raster
            subplot(1,3,2:3)
%             if iROI==4
                PlotRectangle(NbLayers,10,1,0)
%             end
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            subplot(1,3,2:3)
            hold on
            colormap(ColorMap);
            
%             if iROI < 3
%                 Profiles=Profiles*-1;
%             end
            imagesc(mean(imgaussfilt(Profiles,[size(Profiles,1)/100 .0001]),3), CLIM)
%             imagesc(mean(imgaussfilt(Profiles,[0.0001 .0001]),3), CLIM)
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
%             title(NameSwitch{iCdt})
            
            % plot sorting raster
            subplot(1,3,1)
%             if iROI==4
                PlotRectangle(NbLayers,10,1,0)
%             end
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            subplot(1,3,1)
            hold on
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3), CLIM2)
%             imagesc(mean(imgaussfilt(Sorting_Raster,[0.0001 .0001]),3), CLIM2)
            axis([0.5 6.5 0 size(Profiles,1)])
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
%             title('[V-Fix]')
            
            % plot sorting value
            ax = gca;
            YLabel = [];
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, [], [], 0)
            set(gca,'fontsize', 20)
            
            % plot regression coefficient
%             subplot(1,4,4)
%             ax = gca;
%             axis off
%             if iROI < 3
%                 slope=slope*-1;
%             end
%             PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
%                 [0.9 1.6 -0.1 0.5], ToPermute)
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
            
%             SavedTxt = fullfile(FigureFolder,'Attention', ...
%                 ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '.csv']);
%             PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, NameSwitch{iCdt}, slope, ToPermute)
            
            
            fullfile(FigureFolder,'Attention', ...
                ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '.tif'])
            % print
            print(gcf, fullfile(FigureFolder,'Attention', ...
                ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '.tif']), '-dtiff')
            
        end
        
    end
    
    %%
    %     figure('name', ToPlot{iToPlot}, 'Position', [100, 100, 100, 1500], 'Color', [1 1 1], 'Visible', Visibility);
    %     colormap(ColorMap);
    %     imagesc(repmat([CLIM(2):-.01:CLIM(1)]', [1,200]), CLIM)
    %     set(gca,'tickdir', 'out', 'xtick', [],'xticklabel',  [], ...
    %         'ytick', linspace(1,numel(CLIM(2):-.01:CLIM(1)),5),...
    %         'yticklabel', round(10*linspace(CLIM(2),CLIM(1),5))/10, ...
    %         'ticklength', [0.01 0.01], 'fontsize', 8, 'YAxisLocation','right')
    %     box off
    %
    %     print(gcf, fullfile(FigureFolder,'Attention', ...
    %         'GrpLvl_Raster_colorbar.svg'), '-dsvg')
    
end

%% Grp level  ; raster CrossMod = f(Percentile of V stim)
% close all
% mkdir(fullfile(FigureFolder,'CrossSensory'));
% Name = {'AV-A';'AV-V'};
% 
% 
% iRoiGrp = 1;
% 
% for iToPlot = 1:2
%     
%     for iCdt = 1:2
%         
%         for iROI = 1:4
%             
%             if iROI==1
%                 CLIM = [-0.2 0.2];
%                 CLIM2 = CLIM*1.5;
%             else
%                 CLIM = [-0.8 0.8];
%                 CLIM2 = CLIM;
%             end
%             
%             figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
%             
%             NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
%             
%             clear X_sort Profiles Sorting_Raster
%             
%             for iSubj = 1:size(All_X_sort_Cross,1)
%                 Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
%                 X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
%                 Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
%             end
%             
%             X_sort = X_sort(Subj2Include,:);
%             Profiles = Profiles(:,:,Subj2Include);
%             
%             [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
%             
%             
%             % plot main raster
%             subplot(1,3,2:3)
%             PlotRectangle(NbLayers,10,1,0)
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             subplot(1,3,2:3)
%             hold on
%             colormap(ColorMap);
% %             imagesc(mean(imgaussfilt(Profiles,[size(Profiles,1)/200 .0001]),3), CLIM)
%             axis([0.5 6.5 0 size(Profiles,1)])
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             
%             
%             % plot sorting raster
%             subplot(1,3,1)
%             PlotRectangle(NbLayers,10,1,0)
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             subplot(1,3,1)
%             hold on
%             colormap(ColorMap);
% %             imagesc(mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3), CLIM2)
%             axis([0.5 6.5 0 size(Profiles,1)])
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             
%             
%             % plot sorting value
%             ax = gca;
%             YLabel = [];
%             PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, [], [], 0)
%             set(gca,'fontsize', 20)
%             
%             % plot regression coefficient
%             %             subplot(1,4,4)
%             %             ax = gca;
%             %             axis off
%             %             PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
%             %                 [0.9 1.6 -0.1 0.5], ToPermute)
%             %             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
%             
%             
%             SavedTxt = fullfile(FigureFolder,'CrossSensory', ...
%                 ['GrpLvl_Raster_' Name{iCdt} '_fV_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.csv']);
%             PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, [Name{iCdt} '_fV'], slope, ToPermute)
%             
%             
% %             print(gcf, fullfile(FigureFolder,'CrossSensory', ...
% %                 ['GrpLvl_Raster_' Name{iCdt} '_fV_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.tif']), '-dtiff')
%             
%         end
%         
%     end
%     
%     
% end




%% Grp level  ; raster Att = f(Percentile of A)
close all
mkdir(fullfile(FigureFolder,'Attention'));
Name={'Att_A-V', 'Astim_Att_A-V', 'Vstim_Att_A-V'};
NameSwitch={'Att_V-A', 'Astim_Att_V-A', 'Vstim_Att_V-A'};

SubPlot = [3 4 1 2];
iSubPlot = 0;


for iToPlot = 1%:2
    if iToPlot==1
        CLIM = [-.75 .75];
        CLIM2 = CLIM*3;
    else
        CLIM = [-0.4 .4];
        CLIM2 = CLIM*6;
    end
    
    for iCdt = 2 %1:3
        
        for iROI = 1 %1:4
            
            figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_A,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                X_sort(iSubj,:) = All_X_sort_AttfA{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            % plot main raster
%             subplot(1,3,2:3)
%             if iROI==4
%                 PlotRectangle(NbLayers,10,1,0)
%             end
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
            subplot(1,3,2:3)
            hold on
            colormap(ColorMap);
            
            imagesc(mean(imgaussfilt(Profiles,[size(Profiles,1)/100 .0001]),3), CLIM)
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            % plot sorting raster
%             subplot(1,3,1)
% %             if iROI==4
% %                 PlotRectangle(NbLayers,10,1,0)
% %             end
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
            subplot(1,3,1)
            hold on
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3), CLIM2)
            axis([0.5 6.5 0 size(Profiles,1)])
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            % plot sorting value
            ax = gca;
            YLabel = [];
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, [], [], 0)
            set(gca,'fontsize', 20)
            
            % plot regression coefficient
%                         subplot(1,4,4)
%                         ax = gca;
%                         axis off
%                         PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
%                             [0.9 1.6 -0.1 0.5], ToPermute)
%                         set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
            
%             SavedTxt = fullfile(FigureFolder,'Attention', ...
%                 ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '_f(A).csv']);
%             PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, NameSwitch{iCdt}, slope, ToPermute)
%             
            
            
            % print
            print(gcf, fullfile(FigureFolder,'Attention', ...
                ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '.tif']), '-dtiff')
            
        end
        
    end
    
    %%
    %     figure('name', ToPlot{iToPlot}, 'Position', [100, 100, 100, 1500], 'Color', [1 1 1], 'Visible', Visibility);
    %     colormap(ColorMap);
    %     imagesc(repmat([CLIM(2):-.01:CLIM(1)]', [1,200]), CLIM)
    %     set(gca,'tickdir', 'out', 'xtick', [],'xticklabel',  [], ...
    %         'ytick', linspace(1,numel(CLIM(2):-.01:CLIM(1)),5),...
    %         'yticklabel', round(10*linspace(CLIM(2),CLIM(1),5))/10, ...
    %         'ticklength', [0.01 0.01], 'fontsize', 8, 'YAxisLocation','right')
    %     box off
    %
    %     print(gcf, fullfile(FigureFolder,'Attention', ...
    %         'GrpLvl_Raster_colorbar.svg'), '-dsvg')
    
end







%% Grp level  ; raster CrossMod = f(Percentile of V stim)
% close all
% mkdir(fullfile(FigureFolder,'CrossSensory'));
% Name = {'AV-A';'AV-V'};
% 
% 
% iRoiGrp = 1;
% 
% for iToPlot = 1:2
%     
%     for iCdt = 1:2
%         
%         for iROI = 1:4
%             
%             if iROI==1
%                 CLIM = [-0.2 0.2];
%                 CLIM2 = CLIM*1.5;
%             else
%                 CLIM = [-0.8 0.8];
%                 CLIM2 = CLIM;
%             end
%             
%             figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
%             
%             NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
%             
%             clear X_sort Profiles Sorting_Raster
%             
%             for iSubj = 1:size(All_X_sort_Cross,1)
%                 Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,2,iROI};
%                 X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI};
%                 Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
%             end
%             
%             X_sort = X_sort(Subj2Include,:);
%             Profiles = Profiles(:,:,Subj2Include);
%             
%             [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
%             
%             
%             % plot main raster
%             subplot(1,3,2:3)
%             PlotRectangle(NbLayers,10,1,0)
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             subplot(1,3,2:3)
%             hold on
%             colormap(ColorMap);
% %             imagesc(mean(imgaussfilt(Profiles,[size(Profiles,1)/200 .0001]),3), CLIM)
%             axis([0.5 6.5 0 size(Profiles,1)])
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             
%             
%             % plot sorting raster
%             subplot(1,3,1)
%             PlotRectangle(NbLayers,10,1,0)
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             subplot(1,3,1)
%             hold on
%             colormap(ColorMap);
% %             imagesc(mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3), CLIM2)
%             axis([0.5 6.5 0 size(Profiles,1)])
%             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
%                 'ytick', [],'yticklabel', [])
%             
%             
%             % plot sorting value
%             ax = gca;
%             YLabel = [];
%             PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, [], [], 0)
%             set(gca,'fontsize', 20)
%             
%             % plot regression coefficient
%             %             subplot(1,4,4)
%             %             ax = gca;
%             %             axis off
%             %             PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
%             %                 [0.9 1.6 -0.1 0.5], ToPermute)
%             %             set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
%             
%             
%             SavedTxt = fullfile(FigureFolder,'CrossSensory', ...
%                 ['GrpLvl_Raster_' Name{iCdt} '_fA_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.csv']);
%             PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, [Name{iCdt} '_fA'], slope, ToPermute)
%             
%             
% %             print(gcf, fullfile(FigureFolder,'CrossSensory', ...
% %                 ['GrpLvl_Raster_' Name{iCdt} '_fV_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.tif']), '-dtiff')
%             
%         end
%         
%     end
%     
%     
% end