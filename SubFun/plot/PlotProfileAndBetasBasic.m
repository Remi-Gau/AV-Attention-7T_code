function PlotProfileAndBetasBasic(DATA)

if ~exist('DATA', 'var')
    
    Name = 'NoName';
    Data = rand(11,6);
    MVPA = 0;
    Offset = 0;
    
    XY = [0.6 0.2 .04 .07];
    
    %     Betas = randn(11,3);
    %     Thresholds = rand(1,3); %#ok<*NASGU>
    %     WithPerm = 1;
    %     AllNullDist = randn(2048,3);
    
    %     Betas = randn(11,2);
    %     Thresholds = rand(1,2); %#ok<*NASGU>
    %     WithPerm = 1;
    %     AllNullDist = sort(randn(2048,2));
    
    MIN = min(Data(:));
    MAX = max(Data(:));
    
    FigDim = [100 100 1500 1000];
    Visible = 'on';
    
    figure('position', FigDim, 'name', Name, 'Color', [1 1 1], 'visible', Visible)
    
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
    hold on
    grid on
    
else
    
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
    NbSubjects = size(Data,1);
    
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
    
    WithSubj = DATA.WithSubj;
    
    if isfield(DATA, 'MIN')
        MIN = DATA.MIN;
        MAX = DATA.MAX;
    else
        if WithSubj
            MIN = min(Data(:));
            MAX = max(Data(:));
        else
            MIN = min(nanmean(Data)-2*nansem(Data));
            MAX = max(nanmean(Data)+2*nansem(Data));
            if MVPA
                if MIN>0.35
                    MIN = 0.35;
                end
                if MAX<0.65
                    MAX = 0.65;
                end
            else
                if MIN>0
                    MIN = -0.2;
                end
                if MAX<0
                    MAX = 0.2;
                end
            end
        end
    end
    
end

if size(Betas,2)==2;
    WithQuad=0;
else
    WithQuad=1;
end
Scatter = linspace(0,.4,NbSubjects);

if isfield(DATA, 'Color')
    Color = DATA.Color;
else
    Color = 'k';
end

COLOR_Subject= [
    127,127,127;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

COLOR_Layer= [
    254,229,217;
    252,187,161;
    252,146,114;
    251,106,74;
    222,45,38;
    165,15,21];
COLOR_Layer = COLOR_Layer/255;

Transparent = 1;
Fontsize = 13;

NbLayers = size(Data,2);
x = (1:.1:NbLayers)-mean(1:.1:NbLayers);



%%
% box on
hold on
grid on

if MVPA
    plot([0 7], [0.5 0.5], ':k', 'linewidth', 2)
else
    plot([0 7], [0 0], ':k', 'LineWidth', 2)
end

for SubjInd = 1:size(Data,1)
    
    if WithSubj
        plot(1:NbLayers, Data(SubjInd,:), '--', ...
            'LineWidth', 1, 'Color', COLOR_Subject(SubjInd,:), ...
            'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:));
        
        b=Betas(SubjInd,:);
        if WithQuad
            y=b(1)*x+b(2)*x.^2+b(3)+Offset;
        else
            y=b(1)*x+b(2)+Offset; %#ok<*UNRCH>
        end
        %     plot(x-min(x)+1,y,':','LineWidth', 2, 'Color', COLOR_Subject(SubjInd,:))
    end
    
end
if MVPA
    shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
        {'Marker', '.', 'MarkerSize', 25, 'LineStyle', ':', 'LineWidth', 3, 'Color', Color}, Transparent)
else
    shadedErrorBar(1:NbLayers, nanmean(Data),nansem(Data), ...
        {'Marker', '.', 'MarkerSize', 25, 'LineWidth', 3, 'Color', Color}, Transparent)
end

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
        'ticklength', [0.01 0.01], 'fontsize', Fontsize-1)
else
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , 'xticklabel', [], ...
        'ticklength', [0.01 0.01], 'fontsize', Fontsize-1)
end

% t=xlabel('% Cortical depth');
% set(t,'fontsize',Fontsize);

if isfield(DATA,'YLabel')
    t=ylabel(sprintf(DATA.YLabel));
    set(t,'fontsize', Fontsize-1)
end

title(Name,'fontsize', Fontsize)


%% INSETS
XY = DATA.XY;
for i=1:size(Betas,2)
    
    axes('Position',[XY(1)+0.09*(i-1) XY(2) XY(3) XY(4)])
    
    Lim = ceil(max(abs([min(Betas(:,i))-.02 max(Betas(:,i))+.02]))*100)/100;
    
    box off; hold on;
    
    [H(i),P(i)] = ttest(Betas(:,i), 0, 'alpha', 0.05);
    
    for SubjInd=1:size(Betas,1)
        plot(1.2+Scatter(SubjInd), Betas(SubjInd,i), ...
            'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
            'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
    end
    
    plot([0 1.8], [0 0], ':k', 'LineWidth', .5)
    
    h = errorbar(1,nanmean(Betas(:,i)),nansem(Betas(:,i)), '.k');
    
    Sig = [];
    if P(i)<0.001
        Sig = sprintf('ES=%.3f \np<0.001 ',...
            abs(nanmean(Betas(:,i))/nanstd(Betas(:,i))));
    else
        Sig = sprintf('ES=%.3f \np=%.3f ',...
            abs(nanmean(Betas(:,i))/nanstd(Betas(:,i))), P(i));
    end
    
    t = text(.8,Lim+Lim*50/100,sprintf(Sig));
    set(t,'fontsize',Fontsize-2);
    
    if H(i)==1
        set(t,'color','r');
    end
    
    clear Sig
    
    switch i
        case 1
            xTickLabel = 'C';
        case 2
            xTickLabel = 'L';
        case 3
            xTickLabel = 'Q';
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

