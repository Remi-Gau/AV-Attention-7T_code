function PlotInsetFinal(DATA)

ax = DATA.ax;

FontSize = DATA.FontSize;
Betas = DATA.Betas;
ToPermute = DATA.ToPermute;

for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:)';
    Perms(iPerm,:) = mean(Betas.*repmat(tmp2,1,size(Betas,2))); %#ok<AGROW,*SAGROW>
end


for i=1:size(Betas,2)
    
    if isfield(DATA, 'InsetLim')
    Max = DATA.InsetLim(1,i);
    Min = DATA.InsetLim(2,i);
    else
        Max = max(Betas(:));
        Min = min(Betas(:));
    end
    
    
    axes('Position',[...
        ax(1)+(ax(3)/20)+ax(3)/2.15*(i-1) ...
        ax(2) ...
        ax(3)/2.5 ...
        ax(4)*9/10])
    
    box off; hold on;
    
    % plot each subject with a different color
    %         for SubjInd=1:size(tmp,2)
    %             plot(1.2+Scatter(SubjInd), tmp(i,SubjInd), ...
    %                 'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
    %                 'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
    %         end
    
    % plot distribution
    distributionPlot({Betas(:,i)}, 'xValues', 1.3, 'color', [0.8 0.8 0.8], ...
        'distWidth', 0.4, 'showMM', 0, ...
        'globalNorm', 2)
    
    % plot scatter
    h = plotSpread(Betas(:,i), 'distributionIdx', ones(size(Betas(:,i))), ...
        'distributionMarkers',{'o'},'distributionColors',{'w'}, ...
        'xValues', 1.3, 'binWidth', .6, 'spreadWidth', 0.6);
    if ~isnan(h{1})
        set(h{1}, 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...
            'w', 'LineWidth', 1.5)
    end
    
    
    plot([0 1.55], [0 0], ':k', 'LineWidth', .5)

    % plot group average
    errorbar(1,nanmean(Betas(:,i),1),nansem(Betas(:,i),1), 'ok','LineWidth', ...
        1, 'MarkerSize', 1.5);

    
    if isfield(DATA, 'OneSideTTest')
        if strcmp(DATA.OneSideTTest{i},'left')
            P = sum(Perms(:,i)<mean(Betas(:,i)))/numel(Perms(:,i));
        elseif strcmp(DATA.OneSideTTest{i},'right')
            P = sum(Perms(:,i)>mean(Betas(:,i)))/numel(Perms(:,i));
        elseif strcmp(DATA.OneSideTTest{i},'both')
            P = sum(abs(Perms(:,i))>abs(mean(Betas(:,i)))) / numel(Perms(:,i)) ;
        end

    else
        P = sum(abs(Perms(:,i))>abs(mean(Betas(:,i)))) / numel(Perms(:,i)) ;

    end
    
    Sig = []; %#ok<NASGU>
    if P<0.001
        Sig = 'p<0.001';
    else
        Sig = sprintf('p=%.3f',P);
    end
    
    t = text(1.1,Max+(Max-Min)/10,sprintf(Sig));
    set(t,'fontsize',FontSize-1);
    
    if P<0.05
        set(t,'fontsize',FontSize,'fontweight','bold')
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
    
    
    set(gca,'tickdir', 'in', 'xtick', 1.15 , 'xticklabel', xTickLabel, ...
        'ytick', linspace(Min, Max, 5) , 'yticklabel', linspace(Min, Max, 5), ...
        'ticklength', [0.03 0.03], 'fontsize', FontSize-1)
    
    if i==1
        if isfield(DATA,'YLabel')
            t=ylabel(sprintf(DATA.YLabel));
            set(t,'fontsize', FontSize)
        end
    end
    
    if i==2
         set(gca,'yaxislocation', 'right')
    end
    
    axis([0.9 1.55 Min Max])
    
end

end

