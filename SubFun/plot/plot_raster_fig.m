function plot_raster_fig(All_X_sort, All_Profiles, Cdt, Name, ROI_Groups, ROI, MinVert, opt)

iToPlot = opt.iToPlot;
Subj2Include = opt.Subj2Include;
ColorMap = opt.ColorMap;
CLIM = opt.CLIM;
Fontsize = opt.Fontsize;

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = spm_orth(DesMat);


figure('name', opt.name, 'Position', opt.FigDim, 'Color', [1 1 1], 'Visible', opt.Visibility);

for iCdt = 1:size(Cdt, 1)
    
    for iROI = ROI_Groups
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        clear X_sort Profiles Sorting_Raster
        
        for iSubj = 1:size(All_Profiles,1)
            Sorting_Raster(:,:,iSubj) = All_Profiles{iSubj, iToPlot, 1, iROI};
            X_sort(iSubj,:) = All_X_sort{iSubj, iToPlot, iCdt, iROI};
            Profiles(:,:,iSubj) = All_Profiles{iSubj, iToPlot, iCdt, iROI};
        end
        
        X_sort = X_sort(Subj2Include, :);
        Profiles = Profiles(:, :, Subj2Include);
        
        [rho, slope] = CorRegRaster(Profiles, DesMat, iToPlot, X_sort);
        
        subplot( numel(ROI_Groups), size(Cdt, 1), iCdt+size(Cdt, 1)*(iROI-ROI_Groups(1)) )
        hold on
        
        colormap(ColorMap);
        imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), CLIM)
        
        axis([0.5 6.5 0 size(Profiles,1)])
        
        DephLevels = round(linspace(100,0,8));
        DephLevels([1;end]) = [];
        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
            'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
            'ticklength', [0.01 0], 'fontsize', Fontsize)
        
        t=xlabel('cortical depth');
        set(t,'fontsize', Fontsize)
        
        if iROI == ROI_Groups(1)
            title([Name{Cdt(iCdt,2)} ' stimulus'])
        end
        
        ax = gca;
        
        rho = atanh(rho);
        PlotCorrCoeff(ax, slope, opt.name, .15, .1, .05, .05, [0.9 1.4 -1 1], opt.ToPermute);
        
        if iCdt==1
            YLabel = sprintf('%s\nPerc %s %s stim',...
                strrep(ROI(iROI).name,'_','-'), opt.name, Name{Cdt(iCdt,1)});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
        end
        
        clear X
        
    end
    
end

mtit(['Percentile V stim - ' opt.name], 'fontsize', Fontsize, 'xoff', 0, 'yoff', .025)

print(gcf, fullfile(FigureFolder,'Baseline', ...
    ['GrpLvl_Raster_Baseline_V_stim_' ToPlot{iToPlot} '_' ...
    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')

