function GrpBetaCdt = plot_vertex_wise_correlation(BetaCdtion_x, BetaCdtion_y, GrpBetaCdt, NbVertex, VertexWithDataHS, ROI, SubjInd, SubjID, opt)
% outputs group intercept, slope and correlation coefficient of each comparison

Cdt = opt.Cdt;
Name = opt.Name;
ToPlot = opt.ToPlot;
fontsize = opt.fontsize;

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', opt.FigDim, 'Color', [1 1 1], 'Visible', opt.Visibility)
    
    for iCdt = 1:size(Cdt,1)
        
        X_lh = nan(1,NbVertex(1));
        X_lh(1,VertexWithDataHS{1}) = BetaCdtion_x{1}(iToPlot,:,Cdt(iCdt,1));
        X_rh = nan(1,NbVertex(2));
        X_rh(1,VertexWithDataHS{2}) = BetaCdtion_x{2}(iToPlot,:,Cdt(iCdt,1));
        
        Y_lh = nan(1,NbVertex(1));
        Y_lh(1,VertexWithDataHS{1}) = BetaCdtion_y{1}(iToPlot,:,Cdt(iCdt,2));
        Y_rh = nan(1,NbVertex(2));
        Y_rh(1,VertexWithDataHS{2}) = BetaCdtion_y{2}(iToPlot,:,Cdt(iCdt,2));
        
        for iROI = 1:numel(ROI)
            
            X = [X_lh(ROI(iROI).VertOfInt{1}) X_rh(ROI(iROI).VertOfInt{2})];
            Y = [Y_lh(ROI(iROI).VertOfInt{1}) Y_rh(ROI(iROI).VertOfInt{2})];
            
            subplot(numel(ROI), size(Cdt,1), iCdt+size(Cdt,1)*(iROI-1))
            
            hold on
            Beta = PlotScatterDensity(X, Y, opt.Range(iToPlot,:), opt.Range(iToPlot,:), opt.Bins);
            GrpBetaCdt(:,iCdt,iROI,iToPlot,SubjInd) =  Beta;
            
            t=xlabel([Name{Cdt(iCdt,1)} ' stim']);
            set(t,'fontsize', fontsize)
            
            t=ylabel([Name{Cdt(iCdt,2)} ' stim']);
            set(t,'fontsize', fontsize)
            
            t=title(ROI(iROI).name);
            set(t,'fontsize', fontsize)
            
            clear X Y
            
        end
        
        clear Y_lh Y_rh X_lh Y_rh
        
    end
    
    set(gca,'fontsize', fontsize)
    
    mtit(['Subject ' SubjID ' - Baseline - ' ToPlot{iToPlot}], ...
        'fontsize', fontsize, 'xoff', 0, 'yoff', .025)
    
    if opt.print
        print(gcf, fullfile(opt.FigureFolder, ...
            ['Subj_' SubjID '_CorrBaseline_' ToPlot{iToPlot} '.tif']), '-dtiff')
    end
    
end

end