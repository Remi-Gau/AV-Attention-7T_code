function PlotLayersForFig(DATA)


if nargin<1
    error('No data to plot')
    return; %#ok<UNRCH>
end

Fontsize = 12;

Transparent = 1;

% Color for Subjects
COLOR_Subject= [
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130];
COLOR_Subject=COLOR_Subject/255;




%%

Name = strrep(DATA.Name, ' ', '-');
Visible = DATA.Visible;

Mean = DATA.Data.MEAN;
ErrorBar = DATA.Data.SEM;

NbLayers = size(Mean,1);
NbCdts = size(Mean,2);

switch NbCdts
    case 6
        if DATA.MVPA
            m=2; n=3;
        else
            m=3; n=2;
        end
    case 3
        m=1; n=3;
end
SubPlotOrder = DATA.SubPlotOrder;
Legend = DATA.Legend;

if DATA.PlotSub
    Subjects = DATA.Data.grp;
    NbSubjects = size(Subjects,3);
end

if DATA.MVPA
    Beta = DATA.Data.Beta;
else
    Beta = DATA.Data.Beta.DATA;
end

WithQuad = DATA.WithQuad;
Scatter = linspace(0,.4,NbSubjects);

if DATA.PlotSub
    MAX = max(Subjects(:));
    MIN = min(Subjects(:));
else
    MAX = max(Mean(:)+ErrorBar(:));
    MIN = min(Mean(:)+ErrorBar(:));
end

if DATA.MVPA
    XYs = [...
        0.18 0.61 ; ...
        0.465 0.61 ; ...
        0.74 0.61 ; ...
        0.18 0.145 ; ...
        0.465 0.145 ; ...
        0.74 0.145 ; ...
        ];
else
    switch NbCdts
        case 6
            XYs = [...
                0.20 0.72 ; ...
                0.65 0.72 ; ...
                0.20 0.43 ; ...
                0.65 0.43 ; ...
                0.20 0.14 ; ...
                0.65 0.14 ; ...
                ];
        case 3
            XYs = [...
                0.18 0.14 ; ...
                0.465 0.14 ; ...
                0.74 0.14 ; ...
                ];
    end
end


%%
figure('Name', Name, 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);

box off

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for iCdt = 1:NbCdts
    %% Plot main data
    subplot(m,n,iCdt)
    PlotRectangle
    subplot(m,n,iCdt)
    
    hold on; grid on;
    
    shadedErrorBar(1:NbLayers, flipud(Mean(:,SubPlotOrder(iCdt))),flipud(ErrorBar(:,SubPlotOrder(iCdt))), ...
        {'LineWidth', 2, 'Color', 'k'}, Transparent)
    for SubjInd = 1:NbSubjects
        plot(1:NbLayers, flipud(Subjects(:,SubPlotOrder(iCdt),SubjInd)), '--', ...
            'LineWidth', .5, 'Color', COLOR_Subject(SubjInd,:));
    end
    
    if DATA.MVPA
        plot([1 NbLayers], [0.5 0.5], '--k', 'LineWidth', 1)
    else
        plot([1 NbLayers], [0 0], '-k', 'LineWidth', 1)
    end
    
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', linspace(0,100,NbLayers), 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    t=xlabel('Cortical depth');
    set(t,'fontsize',Fontsize);
    
    t=ylabel(Legend{SubPlotOrder(iCdt)}{1});
    set(t,'fontsize',Fontsize+2);
    
    t=title(Legend{SubPlotOrder(iCdt)}{2});
    set(t,'fontsize',Fontsize);
    
    axis([1 NbLayers MIN MAX])
    
    
    %% Inset with betas
    tmp = squeeze(Beta(:,iCdt,:));
    for i=1:size(tmp,1)
        
        Lim = ceil(max(abs([min(tmp(i,:))-.02 max(tmp(i,:))+.02]))*100)/100;
        
        if DATA.MVPA
            XY = XYs(iCdt,:);
            axes('Position',[XY(1)+0.065*(i-1) XY(2) .03 .09])
        else
            switch NbCdts
                case 6
                    XY = XYs(iCdt,:);
                    axes('Position',[XY(1)+0.09*(i-1) XY(2) .05 .07])
                case 3
                    XY = XYs(iCdt,:);
                    axes('Position',[XY(1)+0.075*(i-1) XY(2) .04 .15])
            end
        end
        
        box off; hold on;
        
        [H(i),P(i)] = ttest(tmp(i,:), 0, 'alpha', 0.05);
        
        for SubjInd=1:size(tmp,2)
            plot(1.2+Scatter(SubjInd), tmp(i,SubjInd), ...
                'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
                'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
        end
        
        plot([0 1.8], [0 0], ':k', 'LineWidth', .5)
        
        h = errorbar(1,nanmean(tmp(i,:),2),nansem(tmp(i,:),2), '.k');
        
        Sig = [];
        if P(i)<0.001
            Sig = sprintf('ES=%.3f \np<0.001 ',...
                abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))));
        else
            Sig = sprintf('ES=%.3f \np=%.3f ',...
                abs(nanmean(tmp(i,:))/nanstd(tmp(i,:))), P(i));
        end
        
        t = text(.8,Lim+Lim*50/100,sprintf(Sig));
        set(t,'fontsize',Fontsize-2);
        
        if H(i)==1
            set(t,'color','r');
        end
        
        clear Sig
        
        switch i
            case 1
                xTickLabel = 'L';
            case 2
                if size(tmp,1)==2
                    xTickLabel = 'C';
                else
                    xTickLabel = 'Q';
                end
            case 3
                xTickLabel = 'C';
        end
        
        
        set(gca,'tickdir', 'out', 'xtick', 1.1:3.1 , 'xticklabel', xTickLabel, ...
            'ytick', linspace(Lim*-1, Lim, 5) , 'yticklabel', linspace(Lim*-1, Lim, 5), 'ticklength', [0.05 0.05], 'fontsize', Fontsize-3)
        if i==1
            t=ylabel('betas');
            set(t,'fontsize',Fontsize-3);
        end
        
        axis([0.8 1.8 Lim*-1 Lim])
        
    end
    
    
end

mtit(Name, 'xoff', 0, 'yoff', +0.03, 'fontsize', 16)

print(gcf, fullfile(DATA.FigureFolder, strcat(Name, '_', num2str(NbLayers), 'Layers.pdf')), '-dpdf')
print(gcf, fullfile(DATA.FigureFolder, strcat(Name, '_', num2str(NbLayers), 'Layers.tif')), '-dtiff')

end
