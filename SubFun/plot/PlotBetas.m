function PlotBetas(DATA, Name, Visible, PRINT)

%%
if nargin<1
    error('No data to plot')
    return; %#ok<UNRCH>
end

if nargin<2 || isempty(Name)
    warning('No name specified for this figure')
    Name = 'test';
end

if nargin<3 || isempty(Visible)
    Visible = 'on';
end

if nargin<3 || isempty(PRINT)
    PRINT = 0;
end

Fontsize=9;
% Legend.Bool = 0;
% Legend.Legend = [];

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

COLOR_Subject = COLOR_Subject(DATA.Include,:);

m=3;
n=4;


%% Plot
MAX = DATA.MAX;
MIN = DATA.MIN;

MEAN = DATA.MEAN;

SEM = DATA.SEM;

Legend = DATA.Subjects.Legend;

figure('Name', Name, 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
for CondInd = 1:size(DATA.Subjects.Data,2)
    
    subplot(m,n,CondInd)
    
    hold on
    
    for SubjInd=1:size(DATA.Subjects.Data,3)
        plot(1.2, DATA.Subjects.Data(1,CondInd,SubjInd), ...
            'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
            'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 4)
    end
    
    plot([0 4], [0 0], '--k', 'LineWidth', 2)
    
    errorbar(.8,MEAN(1,CondInd),SEM(1,CondInd), 'MarkerFaceColor',[0 0 1],...
        'MarkerEdgeColor',[0 0 1], 'Marker','o', 'MarkerSize', 4)
    
    t=text(.6, max(DATA.Subjects.Data(1,CondInd,:))+.3, ...
        sprintf('%2.3f', DATA.P(1,CondInd)));
    set(t,'fontsize',Fontsize);
    
    errorbar(1.8,MEAN(2,CondInd),SEM(2,CondInd), 'MarkerFaceColor',[0 0 1],...
        'MarkerEdgeColor',[0 0 1], 'Marker','o', 'MarkerSize', 4)
    for SubjInd=1:size(DATA.Subjects.Data,3)
        plot(2.2, DATA.Subjects.Data(2,CondInd,SubjInd), ...
            'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
            'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 4)
    end
    t=text(1.6, max(DATA.Subjects.Data(2,CondInd,:))+.3, ...
        sprintf('%2.3f', DATA.P(2,CondInd)));
    set(t,'fontsize',Fontsize);
    
    errorbar(2.8,MEAN(3,CondInd),SEM(3,CondInd), 'MarkerFaceColor',[0 0 1],...
        'MarkerEdgeColor',[0 0 1], 'Marker','o', 'MarkerSize', 4)
    for SubjInd=1:size(DATA.Subjects.Data,3)
        plot(3.2, DATA.Subjects.Data(3,CondInd,SubjInd), ...
            'Marker', 'o', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
            'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 4)
    end
    t=text(2.6, max(DATA.Subjects.Data(3,CondInd,:))+.3, ...
        sprintf('%2.3f', DATA.P(3,CondInd)));
    set(t,'fontsize',Fontsize);
    
    axis([0 4 MIN MAX])
    
    set(gca,'tickdir', 'out', 'xtick', 1:3 , ...
        'xticklabel', ['Lin '; 'Quad'; 'Cst '] , 'ticklength', [0.01 0.01], 'fontsize', 10)
    
end

subplot(m,n,1)
t=ylabel('Auditory attention');
set(t,'fontsize',Fontsize);
t=title('Auditory stimulation');
set(t,'fontsize',Fontsize);

subplot(m,n,2)
t=title('Visual stimulation');
set(t,'fontsize',Fontsize);

subplot(m,n,3)
t=title('Audio-Visual stimulation');
set(t,'fontsize',Fontsize);

subplot(m,n,4)
t=title('AV - (A+V)');
set(t,'fontsize',10);

subplot(m,n,5)
t=ylabel('Visual attention');
set(t,'fontsize',Fontsize);

subplot(m,n,9)
if strcmp(Name(1),'V');
    t=ylabel('Visual attention - Auditory attention');
else
    t=ylabel('Auditory attention - Visual attention');
end
set(t,'fontsize',Fontsize);

subplot(m,n,12)
t=text(1.2,MAX-1,['Nb Subject = ' num2str(length(DATA.Include))]);
set(t,'fontsize',Fontsize);

clear Ind t MAX MIN MEAN SEM

if Legend.Bool
    subplot(m,n,12)
    t=legend('Location','BestOutside', Legend.Legend);
    set(t,'fontsize',Fontsize-2);
    legend('boxoff');
end


%%
if PRINT
    switch PRINT
        case 1
            print(gcf, strcat(Name,'.tif'), '-dtiff')
        case 2
            print(gcf, strcat(Name,'.eps'), '-dpsc2')
        case 3
            print(gcf, strcat(Name,'.tif'), '-dtiff')
            print(gcf, strcat(Name,'.eps'), '-dpsc2')
    end
end


end