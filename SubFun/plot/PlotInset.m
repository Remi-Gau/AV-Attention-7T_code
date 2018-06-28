function PlotInset(DATA)

Fontsize = 7;
Betas = DATA.Betas;
ax = DATA.ax;
Thresholds = DATA.Thresholds;
ToPermute = DATA.ToPermute;


for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:)';
    Perms(iPerm,:) = mean(Betas.*repmat(tmp2,1,size(Betas,2))); %#ok<*SAGROW>
end

if size(Betas,2)==2
    
    xPos = [1 1.5];
    
    xLabel = {'Cst','Lin'};
    
    for i=1:size(Betas,2)
        
        if i==1
            ax1 = axes('Position',ax);
        else
            ax2 = axes('Position',ax);
        end
        box off; hold on; grid on;
        
        
        
        if isfield(DATA, 'InsetLim')
            Lim = DATA.InsetLim(i);
        else
            Lim = ceil(max(abs([min(Betas(:,i)) max(Betas(:,i))]))*10)/10;
        end
        
        
%         if isfield(DATA, 'OneSideTTest')
%             P = sum(Perms(x,y,:)>(mean(tmp(x,y,:))-.5))/numel(Perms(x,y,:));
%             [~,P(i)] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i), 'tail', DATA.OneSideTTest{i});
%         else
%             P = sum(Perms(x,y,:)>(mean(tmp(x,y,:))-.5))/numel(Perms(x,y,:));
%             [~,P(i)] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i));
%         end
%         
        if isfield(DATA, 'OneSideTTest')
            if strcmp(DATA.OneSideTTest{i},'left')
                P = sum(Perms(:,i)<mean(Betas(:,i)))/numel(Perms(:,i));
            elseif strcmp(DATA.OneSideTTest{i},'right')
                P = sum(Perms(:,i)>mean(Betas(:,i)))/numel(Perms(:,i));
            elseif strcmp(DATA.OneSideTTest{i},'both')
                P = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(Betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
            end
            
            
%             [~,P] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i), 'tail', DATA.OneSideTTest{i});
            
        else
            
            P = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(Betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
            
            
%             [~,P] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i));
        end
        
        
        
        ToPlot = {Betas(:,i)};
        if isfield(DATA, 'OtherBetas')
            for iOtherBetas = 1:size(DATA.OtherBetas,3)
                ToPlot{end+1} = DATA.OtherBetas(:,i,iOtherBetas);
            end
        end
        
        %     if  i==1
        %         ToPlot = {Betas(:,1)};
        %         for iOtherBetas = 1:size(DATA.OtherBetas,3)
        %             ToPlot{end+1} = DATA.OtherBetas(:,1,iOtherBetas);
        %             %             ToPlot{end+1} = DATA.OtherBetas(:,2,iOtherBetas);
        %         end
        %     else
        %         ToPlot = {Betas(:,2)};
        %         for iOtherBetas = 1:size(DATA.OtherBetas,3)
        %             %             ToPlot{end+1} = DATA.OtherBetas(:,1,iOtherBetas);
        %             ToPlot{end+1} = DATA.OtherBetas(:,2,iOtherBetas);
        %         end
        %     end
        
        distributionPlot(ToPlot, 'xValues', [xPos(i) 3:(1+numel(ToPlot))], 'color', DATA.Color, 'distWidth', 0.45, 'showMM', 0, ...
            'globalNorm', 2)
        
        h = plotSpread(Betas(:,i), 'distributionIdx', ones(size(Betas(:,i))), ...
            'distributionMarkers',{'o'},'distributionColors',{'w'}, ...
            'xValues', xPos(i), 'binWidth', 0.2, 'spreadWidth', 0.4);
        set(h{1}, 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        
        plot([0 3], [0 0], ':k', 'LineWidth', 2)
        
        Sig = '';
        if P<.01
            Sig = [Sig ' **'];
        elseif P<.05
            Sig = [Sig ' *'];
        else
        end
        
        xLabel{i} = [xLabel{i} Sig];
        
        axis([0.7 1.8 Lim*-1-Lim*20/100 Lim+Lim*20/100])
        
        set(gca,'tickdir', 'out', ...
            'ytick', linspace(Lim*-1, Lim, 5) , 'yticklabel', {Lim*-1, '', 0, '', Lim}, ...
            'ticklength', [0.02 0.02], 'fontsize', Fontsize-2,...
            'XGrid', 'off')
        
        if i==2
            set(ax2,'Color','none');
            set(ax2,'Xtick',[]);
            set(ax2,'YAxisLocation','right',...
                'XGrid', 'off')
            
            set(ax1,'xtick', xPos , 'xticklabel', xLabel,...
                'XGrid', 'off')
            
        else
            if DATA.YLabelInset
                t=ylabel(sprintf('Param. est. [a u]'));
                set(t,'fontsize', Fontsize-1)
            end
        end
        
    end
    
    
    
    
elseif size(Betas,2)==3
    
    for i=1:size(Betas,2)
        
        Lim = ceil(max(abs(Betas(:,i)))*10)/8;
        
        
        axes('Position',[ax(1)+(0.01+ax(3)/3)*(i-1) ax(2) ax(3)/4 ax(4)])
        
        box off; hold on;
        

        
        %         for SubjInd=1:size(tmp,2)
        %             plot(1.2+Scatter(SubjInd), tmp(i,SubjInd), ...
        %                 'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
        %                 'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
        %         end
        
        distributionPlot({Betas(:,i)}, 'xValues', 1.3, 'color', [0.8 0.8 0.8], ...
            'distWidth', 0.4, 'showMM', 0, ...
            'globalNorm', 2)
        
        h = plotSpread(Betas(:,i), 'distributionIdx', ones(size(Betas(:,i))), ...
            'distributionMarkers',{'o'},'distributionColors',{'w'}, ...
            'xValues', 1.3, 'binWidth', .5, 'spreadWidth', 0.5);
        if ~isnan(h{1})
            set(h{1}, 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1)
        end
        
        plot([0 1.55], [0 0], ':k', 'LineWidth', .5)
        
        h = errorbar(1,nanmean(Betas(:,i),1),nansem(Betas(:,i),1), '.k');
        
        
        
        if isfield(DATA, 'OneSideTTest')
            if strcmp(DATA.OneSideTTest{i},'left')
                P = sum(Perms(:,i)<mean(Betas(:,i)))/numel(Perms(:,i));
            elseif strcmp(DATA.OneSideTTest{i},'right')
                P = sum(Perms(:,i)>mean(Betas(:,i)))/numel(Perms(:,i));
            elseif strcmp(DATA.OneSideTTest{i},'both')
                P = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(Betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
            end


%             [~,P] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i), 'tail', DATA.OneSideTTest{i});
            
        else
            P = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(Betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;     

%             [~,P] = ttest(Betas(:,i), 0, 'alpha', Thresholds(i));

        end
        
        Sig = [];
        if P<0.001
            Sig = 'p<0.001';
        else
            Sig = sprintf('p=%.3f',P);
        end
                
        t = text(1,Lim,sprintf(Sig));
        set(t,'fontsize',Fontsize-1);
        
        if P<0.05
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
        
        
        set(gca,'tickdir', 'in', 'xtick', 1.3 , 'xticklabel', xTickLabel, ...
            'ytick', linspace(Lim*-1, Lim, 5) , 'yticklabel', linspace(Lim*-1, Lim, 5), ...
            'ticklength', [0.03 0.03], 'fontsize', Fontsize-2)
        
        if i==1
            t=ylabel(sprintf('Param. est. [a u]'));
            set(t,'fontsize', Fontsize)
        end
        
        axis([0.9 1.55 Lim*-1 Lim])
        
    end
    
end

end

