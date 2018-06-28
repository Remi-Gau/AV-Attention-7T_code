function PlotCorrCoeff(ax, Values2Plot, XTickLabel, XOffSet, YOffSet, XDim, YDim, Axis, ToPermute)

if nargin<9
    ToPermute=[];
end

axPos = ax.Position;
axPos(1) = axPos(1)+XOffSet;
axPos(2) = axPos(2)+YOffSet;
axPos(3) = XDim;
axPos(4) = YDim;
axes('Position',axPos);
hold on



if isempty(ToPermute)
    [~,P] = ttest(Values2Plot);
else
    Perms = ToPermute.*repmat(Values2Plot,[size(ToPermute,1),1]);
    Perms = mean(Perms,2);
    P = sum( abs(Perms) > abs( mean(Values2Plot) ) ) / numel(Perms);
end


h = errorbar(1, mean(Values2Plot), nansem(Values2Plot), 'ok','LineStyle','none','Color',[0 0 0]);
set(h, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5)

distributionPlot(Values2Plot', 'xValues', 1.3, 'color', [0.8 0.8 0.8], ...
    'distWidth', 0.4, 'showMM', 0, ...
    'globalNorm', 2)


h = plotSpread(Values2Plot', ...
    'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
        'distributionMarkers',{'o'},'distributionColors',{'w'}, ...
        'xValues', 1.3, 'binWidth', .5, 'spreadWidth', 0.5);
set(h{1}, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5)


plot([0.9 1.9],[0 0],':k', 'linewidth', 2)

ax = axis;
MAX = max(abs([ax(3) ax(4)]));
Y_scale = linspace(MAX*-1,MAX,5);
set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
    'xtick', 1.15 ,'xticklabel',XTickLabel, 'ytick', Y_scale ,'yticklabel', Y_scale,...
    'ygrid', 'on','fontsize',10)

axis([Axis(1) Axis(2) ax(3) ax(4)])


if P<0.001
    Sig = sprintf('\np<0.001 ');
else
    Sig = sprintf('\np=%.3f ',P);
end
t = text(1.1,ax(4),sprintf(Sig));
set(t,'fontsize',10);

if P<.05
%     set(t,'color','r');
end

% t=ylabel('Param. est. [a u]');
% set(t,'fontsize',10)

end