function GridSearchPlot(grid, name, opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
       
    figure('name', name, 'Color', [1 1 1]);
    
    subplot(211)
    
    box off
    hold on
    
    plot(grid.first.acc, 'linewidth', 2)
    
    axis([0.5 size(grid.first.acc,2)+.5 ...
        min(grid.first.acc)-10 max(grid.first.acc)+10])
    
    set(gca,'tickdir', 'out', ...
        'xtick', 1:size(grid.first.acc,2), ...
        'xticklabel',  0:size(grid.first.acc,2)-1, ...
        'ytick', 0:10:100, ...
        'yticklabel', 0:10:100)
    
     title('1^{rst} grid search');
     
     xlabel('c parameters values (powers of 2)');
     ylabel('Percent correct');
     
     [val, id] = max(grid.first.acc);
     
     plot(id, val, 'or', 'MarkerSize', 8, 'linewidth', 2)
     
     plot([id-1 id+1], [val+9 val+9], 'k', 'linewidth', 3)
     
    subplot(212)
    
    box off
    hold on
    
    plot(grid.second.acc, 'linewidth', 2)
    
    axis([0.5 size(grid.second.acc,2)+.5 ...
        min(grid.second.acc)-10 max(grid.second.acc)+10])
    
    set(gca,'tickdir', 'out', ...
        'xtick', 1:size(grid.second.acc,2), ...
        'xticklabel',  grid.second.param, ...
        'ytick', 0:10:100, ...
        'yticklabel', 0:10:100)
    
     title('2^{nd} grid search');
     
     xlabel('c parameters values');
     ylabel('Percent correct');
     
     [val, id] = max(grid.second.acc);
     
     plot(id, val, 'or', 'MarkerSize', 8, 'linewidth', 2)
     
     if opt.print.do
          print(gcf, fullfile(opt.print.folder, [name '.tif']), '-dtiff')
     end


end

