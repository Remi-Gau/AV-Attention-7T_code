function plot_group_vertex_correlation(GrpBetaCdt, ROI, CdtionName, opt)

Xpos = opt.Xpos;
fontsize = opt.fontsize;

figure('name', opt.name, 'Position', opt.FigDim, 'Color', [1 1 1], 'Visible', opt.Visibility)

for iCdt=1:numel(CdtionName)
    
    for iROI=1:1:numel(ROI)
        
        subplot(numel(ROI), numel(CdtionName), iCdt+numel(CdtionName)*(iROI-1))
        hold on
        
        plot([0.9 1.9],[0 0],':k', 'linewidth', 2)
        
        tmp = squeeze(GrpBetaCdt(3, iCdt, iROI, :, opt.Subj2Include));
        tmp = atanh(tmp);
        for i=1:size(tmp,1)
            [P,H] = SignPermTest(tmp(i,:)');
            if P<0.001
                Sig = sprintf('\np<0.001 ');
            else
                Sig = sprintf('\np=%.3f ',P);
            end
            t = text(Xpos(i), 1.75, sprintf(Sig));
            set(t,'fontsize', fontsize);
            
            if H==1
                set(t,'color','r');
            end
        end
        h = errorbar(Xpos, mean(tmp,2), nanstd(tmp,2), 'o','LineStyle','none','Color',[0 0 0]);
        set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        h = plotSpread(tmp', ...
            'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
            'xValues', [1.2 1.7], 'binWidth', 0.25, 'spreadWidth', 0.6);
        set(h{1}, 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1)
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', fontsize, ...
            'xtick', [1.1 1.6] ,'xticklabel', {'Cst','Lin'}, 'ytick', -3:.25:3 ,'yticklabel', -3:.25:3,...
            'ygrid', 'on')
        axis([0.9 1.8 -0.5 1.75])
    end
end

for iCdt=1:numel(CdtionName)
    subplot(numel(ROI),numel(CdtionName), iCdt)
    t=title(CdtionName{iCdt});
    set(t,'fontsize',fontsize)
end

for iROI=1:1:numel(ROI)
    subplot(numel(ROI), numel(CdtionName), 1+numel(CdtionName)*(iROI-1))
    t=ylabel(sprintf('%s \n Corr Coef\n(Fishcer trans)',ROI(iROI).name));
    set(t,'fontsize', fontsize)
end

if opt.print
    print(gcf, fullfile(opt.FigureFolder, [opt.name '.tif']), '-dtiff')
end

end