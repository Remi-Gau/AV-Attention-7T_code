function PlotProfileAndBetas(DATA)

if ~exist('DATA', 'var')
    
    DATA.Name = 'NoName';
    
    DATA.MVPA = 0;
    DATA.Data = rand(11,6);
    DATA.Betas = randn(11,2);
    
    DATA.Thresholds = rand(1,2); %#ok<*NASGU>
    
    DATA.WithSubj = 1;
    DATA.WithPerm = 1;
    
    FigDim = [100 100 1500 1000];
    Visible = 'on';
    
    figure('position', FigDim, 'name', DATA.Name, 'Color', [1 1 1], 'visible', Visible)
    
end

if isfield(DATA, 'Name')
    Name = DATA.Name;
else
    Name = 'NoName';
end

if isfield(DATA, 'Data')
    Data = DATA.Data;
else
    error('Got no subject profile data to plot')
end

if isfield(DATA, 'Betas')
    Betas = DATA.Betas;
else
    error('Got no subject beta data to plot')
end

MVPA = DATA.MVPA;
if MVPA
    Offset = .5;
else
    Offset = 0;
end

if isfield(DATA, 'Color')
    Color = DATA.Color;
else
    Color = 'k';
end

if isfield(DATA, 'MIN')
    MIN = DATA.MIN;
    MAX = DATA.MAX;
else
    MIN = min(Data(:));
    MAX = max(Data(:));
end

if isfield(DATA, 'Transparent')
    Transparent = DATA.Transparent;
else
    Transparent = 1;
end

if isfield(DATA, 'FontSize')
    FontSize = DATA.FontSize;
else
    FontSize = 12;
end

WithSubj = DATA.WithSubj;


if size(Betas,2)==2
    WithQuad=0;
else
    WithQuad=1;
end

COLOR_Subject = ColorSubject();

COLOR_Layer= [
    254,229,217;
    252,187,161;
    252,146,114;
    251,106,74;
    222,45,38;
    165,15,21];
COLOR_Layer = COLOR_Layer/255;


NbLayers = size(Data,2);
x = (1:.1:NbLayers)-mean(1:.1:NbLayers);



%%
hold on
grid on

% plot zero line
if MVPA
    plot([0 10], [0.5 0.5], ':k', 'linewidth', 2)
else
    plot([0 10], [0 0], ':k', 'LineWidth', 2)
end

for SubjInd = 1:size(Data,1)
    
    if WithSubj
        
        % to plot each particpipant with a different color
        %         plot(1:NbLayers, Data(SubjInd,:), '-', ...
        %             'LineWidth', 1.5, 'Color', COLOR_Subject(SubjInd,:), ...
        %             'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:));
        
        plot(1:NbLayers, Data(SubjInd,:), '-', ...
            'LineWidth', 1, 'Color', [.8 .8 .8], ...
            'Marker', '.', 'MarkerEdgeColor', [.8 .8 .8]);
        
        % computes predicted profile of each subject
        b=Betas(SubjInd,:);
        if WithQuad
            y=b(1)+b(2)*x+b(3)*x.^2+Offset;
        else
            y=b(1)+b(2)*x+Offset; %#ok<*UNRCH>
        end
        
        % plot predicted of each subject
        %         plot(x-min(x)+1,y,'-', ...
        %             'LineWidth', 1.5, 'Color', COLOR_Subject(SubjInd,:));
        
        %         plot(x-min(x)+1,y,'-', ...
        %             'LineWidth', 1, 'Color', [.7 .7 .7]);
    end
    
end

% plot the group mean
if MVPA
    shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
        {'Marker', 'x', 'MarkerSize', 8, 'LineStyle', '--', 'LineWidth', 2, 'Color', Color}, Transparent)
else
    shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
        {'Marker', '.', 'MarkerSize', 25, 'LineWidth', 3, 'Color', Color}, Transparent)
end

% to plot the predicted qt the group level
b = mean(Betas);
if WithQuad
    y = x*(b(1))+b(2)*x.^2+b(3)+Offset;
else
    y = x*(b(1))+b(2)+Offset;
end
% plot(x-min(x)+1,y,':k','LineWidth', 2)

axis([0.5 NbLayers+.5 MIN MAX])

if MVPA
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , 'ytick', 0:.1:1, ...
        'yticklabel', 100*0:.1:1, 'xticklabel', [], ...
        'ticklength', [0.01 0.01], 'fontsize', FontSize-1,...
        'XGrid', 'off')
else
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , 'xticklabel', [], ...
        'ticklength', [0.01 0.01], 'fontsize', FontSize-1,...
        'XGrid', 'off')
end

if isfield(DATA,'YLabel')
    t=ylabel(sprintf(DATA.YLabel));
    set(t,'fontsize', FontSize-1)
end

title(Name,'fontsize', FontSize)

ax=axis;

end

