function plot_group_vertex_slope_corr(GrpBetaCdt, CdtionName, Cdtion, ROI, opt)

Xpos = opt.Xpos;
fontsize = opt.fontsize;

SubPlot = [1 3 2 4];

COLOR_Subject = ColorSubject();

for iToPlot=2:3
    
    figure('name', opt.name, 'Position', opt.FigDim, 'Color', [1 1 1], 'Visible', opt.Visibility)
    
    switch iToPlot
        case 2
            suffix='Slope';
            YLabel = 'Slope';
        case 3
            suffix='CorrCoeff';
            YLabel = 'Corr Coef\n(Fishcer trans)';
    end
    
    for iROI=1:1:numel(ROI)
        
        if iROI==1
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        elseif iROI==3
            SubPlotTitle = CdtionName{Cdtion(iROI)};
        else
            SubPlotTitle = [];
        end
        
        subplot(2, 2, SubPlot(iROI))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCdt(iToPlot, Cdtion(iROI), iROI, :, opt.Subj2Include));
        if iToPlot==3
            tmp = atanh(tmp);
        end
        
        for i=1:size(tmp,1)
            [P,H] = SignPermTest(tmp(i,:)');
            if P<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P);
            end
            t = text(Xpos(i), .5, sprintf(Sig));
            set(t,'fontsize', fontsize);
            
            if H==1
                set(t,'color','r');
            end
        end

        
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        for i=1:size(tmp,2)
            plot(1.025+i*.025,tmp(1,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
            plot(1.525+i*.025,tmp(2,i),'.','color', COLOR_Subject(i,:), 'MarkerSize', 14)
        end
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', fontsize, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.5 0.5])
        
        t=title(sprintf('%s\n%s', SubPlotTitle, ROI(iROI).name));
        set(t,'fontsize', fontsize)
        
        t=ylabel(sprintf(YLabel));
        set(t,'fontsize', fontsize)
    end
    
    mtit(strrep([opt.name '-' suffix], '_', '-'), 'fontsize', fontsize, 'xoff',0,'yoff',.05)
    
    if opt.print
        print(gcf, fullfile(FigureFolder, [opt.name '_' suffix '.tif']), '-dtiff')
    end
    
end



end