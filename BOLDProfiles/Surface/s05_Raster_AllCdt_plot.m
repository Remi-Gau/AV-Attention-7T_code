
%%
% Color map
X = 0:0.001:1;
R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
ColorMap = [R' G' B'];
clear R G B

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

StartFolder='/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

load(fullfile(StartFolder,'Results','Profiles','Surfaces','RasterAllCdt_CV.mat'))

Subj2Include = true(11,1);

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

ToPlot={'Cst','Lin'};

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')


CLIM = [-5 5];


%% Grp level  ; raster stim = f(Percentile of V)
close all
Cdt =[2 1;2 2;2 3;2 4;2 5;2 6];
Name={'A','V','AV','A','V','AV'};
ToPlot={'Cst','Lin'};

for iToPlot = 1%:numel(ToPlot)
    
    close all
    
    for iROI = 1:4 %:numel(ROI)
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        fprintf('    %s\n',ROI(iROI).name)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
        
        for iCdt = 1:6
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_V,1)
                Sorting_Raster(:,:,iSubj) = (All_Profiles_V{iSubj,iToPlot,2,iROI}+All_Profiles_V{iSubj,iToPlot,5,iROI})/2;
                X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            MeanProfiles = mean(Profiles,3);
            
            
            subplot(2,3,iCdt)
            
            colormap(ColorMap);
            imagesc(flipud(imgaussfilt(MeanProfiles,[10 .001])), CLIM)
            %             imagesc(MeanProfiles, [-5 5])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t = title(strrep(Conditions_Names{iCdt},'ention', ''));
            set(t,'fontsize',10)
            
            
            ax = gca;
            axes('Position',ax.Position);
            hold on
            errorbar(1:6,mean(MeanProfiles),nansem(squeeze(mean(Profiles)),2), 'k', 'linewidth', 2)
            for iSubj=1:size(Profiles,3)
                plot(1:6,mean(Profiles(:,:,iSubj)), ':k', 'linewidth', .5)
            end
            plot([1 6],[0 0], '--k')
            
            axis([0.5 6.5 -3 6])
            
            DephLevels = round(linspace(100,0,8));
            DephLevels([1;end]) = [];
            set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                'YAxisLocation','right', 'ytick', -10:10,'yticklabel', -10:10, ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            
            if iCdt==1 ||  iCdt==4
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s %s stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
            end
            
            if iCdt==3 ||  iCdt==6
                PlotColorBar(ax, ColorMap, CLIM)
            end
        end
        
        mtit([strrep(ROI(iROI).name,'_',' ') ' - Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_RasterAllCdt_Baseline_V_stim_' ToPlot{iToPlot} '_' ROI(iROI).name '_CV.tiff']), '-dtiff')
        
    end
    
end


%% Grp level  ; raster stim = f(Percentile of A)
close all
Cdt =[1 1;1 2;1 3;1 4;1 5;1 6];
Name={'A','V','AV','A','V','AV'};
ToPlot={'Cst','Lin'};

for iToPlot = 1%:numel(ToPlot)
    
    close all
    
    for iROI = 1:4 %numel(ROI)
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        fprintf('    %s\n',ROI(iROI).name)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
        
        for iCdt = 1:6
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_V,1)
                Sorting_Raster(:,:,iSubj) = (All_Profiles_A{iSubj,iToPlot,1,iROI}+All_Profiles_V{iSubj,iToPlot,4,iROI})/2;
                X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
            end
            
            Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            MeanProfiles = mean(Profiles,3);
            
            
            subplot(2,3,iCdt)
            
            colormap(ColorMap);
            imagesc(flipud(imgaussfilt(MeanProfiles,[10 .001])), CLIM)
            %             imagesc(MeanProfiles, [-5 5])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t = title(strrep(Conditions_Names{iCdt},'ention', ''));
            set(t,'fontsize',10)
            
            
            ax = gca;
            axes('Position',ax.Position);
            hold on
            errorbar(1:6,mean(MeanProfiles),nansem(squeeze(mean(Profiles)),2), 'k', 'linewidth', 2)
            for iSubj=1:size(Profiles,3)
                plot(1:6,mean(Profiles(:,:,iSubj)), ':k', 'linewidth', .5)
            end
            plot([1 6],[0 0], '--k')
            
            axis([0.5 6.5 -3 6])
            
            DephLevels = round(linspace(100,0,8));
            DephLevels([1;end]) = [];
            set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                'YAxisLocation','right', 'ytick', -10:10,'yticklabel', -10:10, ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            
            if iCdt==1 ||  iCdt==4
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s %s stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
            end
            
            if iCdt==3 ||  iCdt==6
                PlotColorBar(ax, ColorMap, CLIM)
            end
            
        end
        
        mtit([strrep(ROI(iROI).name,'_',' ') ' - Percentile A stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_RasterAllCdt_Baseline_A_stim_' ToPlot{iToPlot} '_' ROI(iROI).name '_CV.tiff']), '-dtiff')
        
    end
    
end
