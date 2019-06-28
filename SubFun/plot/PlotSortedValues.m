function PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, PlotScale, plot_sub_zero)


if nargin<7
    plot_sub_zero=1;
end

MAX = ceil(max(abs(mean(X_sort,1))));

axPos = ax.Position;
axPos(1) = axPos(1)-.06;
axPos(3) = .05;

% plot the sorting variable with horizontal error bars
axes('Position',axPos);
hold on

hh = herrorbar(mean(X_sort,1),1:NbBin, nanstd(X_sort,1));
set(hh,  'color', 'k', 'linewidth',.5)
plot([0 0],[0 NbBin],'-k', 'linewidth',.5)

% plot each subject's zero
if plot_sub_zero
    for i=1:size(X_sort,1)
        tmp = abs(X_sort(i,:));
        [~, idx] = min(tmp); %index of closest value
        plot([MAX*-1 MAX],[idx idx],':k','linewidth',2)
    end
end

axis([MAX*-1 MAX 1 size(Profiles,1)])

XScale = mat2cell([MAX*-1 0 MAX], 1, [1 1 1]);
XScale = cellfun(@num2str,XScale,'UniformOutput',0);
if PlotScale
    set(gca,'color', 'none', 'tickdir', 'out', 'xtick', [MAX*-1 0 MAX],'xticklabel', XScale, ...
        'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
        'ticklength', [0.01 0.01], 'fontsize', 10)
        set(gca,'color', 'none', 'tickdir', 'out', 'xtick', [MAX*-1 0 MAX],'xticklabel', XScale, ...
        'ytick', linspace(1,size(Profiles,1),5),'yticklabel', 0:25:100, ...
        'ticklength', [0.01 0.01], 'fontsize', 10)
else
    set(gca,'color', 'none', 'tickdir', 'out', 'xtick', [MAX*-1 0 MAX],'xticklabel', [MAX*-1 0 MAX], ...
    'ytick', linspace(1,size(Profiles,1),5),'yticklabel', [], ...
    'ticklength', [0.01 0.01], 'fontsize', 10)
end

t=ylabel(YLabel);
set(t,'fontsize',10)

end
