function PlotROIForFig(DATA)

% Color for Subjects
COLOR_Subject= [
    0,0,0;
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

Fontsize = 12;

Mean = nanmean(DATA.Data);
ErrorBar = nansem(DATA.Data);

ToPermute = DATA.ToPermute;

for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:)';
    Perms(iPerm,:) = mean(DATA.Data.*repmat(tmp2,1,size(DATA.Data,2))); %#ok<*SAGROW>
end

MAX = DATA.MAX;
MIN = DATA.MIN;

Legend = DATA.Legend;

NbSubjects = numel(DATA.Data);
Scatter = linspace(0,0.1,NbSubjects);

hold on; grid on;

errorbar(1, Mean,ErrorBar, 'o', 'LineWidth', 1, 'Color', 'k')

% h = plotSpread(DATA.Data, 'distributionIdx', ones(size(DATA.Data)), ...
%     'distributionMarkers',{'o'},'distributionColors',{'k'}, ...
%     'xValues', 1.2, 'binWidth', 0.1, 'spreadWidth', 0.2);
% set(h{1}, 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1)

for SubjInd = 1:NbSubjects
    plot(1.05+Scatter(SubjInd), DATA.Data(SubjInd), ...
        'linestyle', 'none', ...
        'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
        'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 28)
end

if DATA.MVPA
    plot([0 2], [0.5 0.5], '--k', 'LineWidth', 1)
else
    plot([0 2], [0 0], '-k', 'LineWidth', 1)
end

if DATA.MVPA
    [H,P] = ttest(DATA.Data-.5, 0, 'alpha', 0.05);
    ES =  abs(nanmean(DATA.Data-.5)/nanstd(DATA.Data-.5));
else
    
%     if isfield(DATA, 'OneSideTTest')
%         [H,P] = ttest(DATA.Data, 0, 'alpha', 0.05, 'tail', DATA.OneSideTTest);
%     else
%     [H,P] = ttest(DATA.Data, 0, 'alpha', 0.05);
%     end
    
    if isfield(DATA, 'OneSideTTest')
        if strcmp(DATA.OneSideTTest,'left')
            P = sum(Perms<mean(DATA.Data))/numel(Perms);
        elseif strcmp(DATA.OneSideTTest,'right')
            P = sum(Perms>mean(DATA.Data))/numel(Perms);
        elseif strcmp(DATA.OneSideTTest,'both')
            P = sum( abs((Perms-mean(Perms))) > abs((mean(DATA.Data)-mean(Perms))) ) / numel(Perms) ;
        end
    else
        P = sum( abs((Perms-mean(Perms))) > abs((mean(DATA.Data)-mean(Perms))) ) / numel(Perms) ;
    end
    
    
    ES =  abs(nanmean(DATA.Data)/nanstd(DATA.Data));
end


Sig = [];
if P<0.001
    Sig = sprintf('ES=%.3f \np<0.001 ',ES);
else
    Sig = sprintf('ES=%.3f \np=%.3f ',ES, P);
end

t = text(1,MAX-MAX*25/100,sprintf(Sig));
set(t,'fontsize',Fontsize-2);

if P<.05
    set(t,'color','r');
end

clear Sig


set(gca,'tickdir', 'out', 'xtick', [] , ...
    'xticklabel', [], 'ticklength', [0.01 0.01], 'fontsize', Fontsize)

t=ylabel(Legend{1});
set(t,'fontsize',Fontsize+2);

t=title(Legend{2});
set(t,'fontsize',Fontsize);

axis([0.95 1.2 MIN MAX])



end
