function PlotLayers(DATA, Name, Visible, PRINT)

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

Fontsize = 12;
Transparent = 0;
Legend.Bool = 0;
Legend.Legend = [];

% Color for Conditions
COLOR =   [...
    255 150 150; ...
    150 255 150; ...
    150 220 255; ...
    150 150 150; ...
    255 75 75; ...
    75 255 75; ...
    75 75 255; ...
    75 75 75; ...
    255 0 0; ...
    0 255 0; ...
    0 0 255; ...
    0 0 0];
COLOR = repmat([0 0 0], size(COLOR,1),1);
COLOR=COLOR/255;

% Color for Subjects
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

% Color for sub-ROI
COLOR_SubROI = [
    0 0 1; ...
    1 0 0; ...
    0 1 0; ...
    1 0 1; ...
    1 1 0; ...
    0.5 .5 0.5; ...
    0 0 0];

% Color for blocks
COLOR_Blocks= [
    0,0,250;
    0,0,250;
    0,0,250;
    250,0,0;
    250,0,0;
    250,0,0;
    0,250,0;
    0,250,0;
    0,250,0;
    0,0,0;
    0,0,0;
    0,0,0;    ];
COLOR_Blocks=COLOR_Blocks/255;


%%
MAX = DATA.MAX;

MIN = DATA.MIN;

Name = strrep(Name, ' ', '-');

if isfield(DATA, 'shadedErrorBar')
    Mean = DATA.shadedErrorBar.mean;
    ErrorBar = DATA.shadedErrorBar.errorbar;
    NbLayers = size(Mean,1);
    NbConditions = size(Mean,2);
    try
        NbSubROIs = size(Mean,3);
        Legend = DATA.shadedErrorBar.Legend;
    catch
    end
end

if isfield(DATA, 'Blocks')
    BlocksMeans = DATA.Blocks.mean;
    BlocksErrorBar = DATA.Blocks.errorbar;
end

if isfield(DATA, 'Subjects')
    Subjects = DATA.Subjects.Data;
    NbSubjects = size(Subjects,3);
    SubjectList = DATA.Subjects.List;
    Transparent = 1;
    Legend = DATA.Subjects.Legend;
    COLOR_Subject = COLOR_Subject(DATA.Include,:);
    SubPlots = DATA.Subjects.SubPlots;
end

if isfield(DATA, 'Min')
    Min = DATA.Min;
    Max = DATA.Max;
end

if isfield(DATA, 'P')
    P = DATA.P;
end

if isfield(DATA, 'BOLD_Thresholds')
    BOLD_Thresholds = DATA.BOLD_Thresholds;
end

if NbConditions == 6
    m=2; n=3;
    COLOR = COLOR([1 2 3 5 6 7],:);
elseif NbConditions == 2
    m=1; n=2;
    COLOR =   [...
        180 0 180; ...
        0 180 180];
    COLOR=COLOR/255;
elseif NbConditions == 4
    m=2; n=2;
    COLOR =   [...
        75 255 75; ...
        75 75 255; ...
        180 0 180; ...
        0 180 180];
    COLOR=COLOR/255;
elseif NbConditions == 10
    COLOR =   [...
        255 150 150; ...
        150 255 150; ...
        150 220 255; ...
        70 70 105; ...
        70 35 105; ...
        255 75 75; ...
        75 255 75; ...
        75 75 255; ...
        180 0 180; ...
        0 180 180];
    COLOR=COLOR/255;
    
    m=2; n=5;
else
    m=3; n=4;
end

if isfield(DATA, 'MainEffect')
    MainEffect = DATA.MainEffect;
else
    MainEffect = 0;
end

%%
figure('Name', Name, 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for CondInd=1:NbConditions
    
    subplot(m,n,CondInd)
    
    hold on; grid on; box on;
    
    if exist('Subjects', 'var')
        if any(CondInd==SubPlots)
            for SubjInd = 1:NbSubjects
                plot(1:NbLayers, flipud(Subjects(:,CondInd,SubjInd)), '--', ...
                    'LineWidth', 1, 'Color', COLOR_Subject(SubjInd,:), ...
                    'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:));
            end
            clear SubjInd
        end
    end
    
    if exist('NbSubROIs', 'var')
        if NbSubROIs>1
            for SubROI_Ind = 1:NbSubROIs
                errorbar([1:NbLayers]+.1*SubROI_Ind, ...
                    Mean(:,CondInd,SubROI_Ind), ...
                    ErrorBar(:,CondInd,SubROI_Ind), ...
                    'Color', COLOR_SubROI(SubROI_Ind,:),...
                    'LineWidth', 1)
            end
            
        else
            shadedErrorBar(1:NbLayers, Mean(:,CondInd),ErrorBar(:,CondInd), ...
                {'Marker', '.', 'MarkerSize', 28, 'LineWidth', 3, 'Color', COLOR(CondInd,:)}, Transparent)
        end
    end
    
    
    plot([0.5 NbLayers+.5], [0 0], '-k', 'LineWidth', 2)
    
    if exist('BlocksMeans', 'var')
        for BlockInd = 1:size(BlocksMeans,2)
            errorbar(...
                [1:NbLayers], ...
                flipud(BlocksMeans(:,BlockInd,CondInd)), ...
                flipud(BlocksErrorBar(:,BlockInd,CondInd)), ...
                'Color', COLOR_Blocks(BlockInd,:)) %#ok<NBRAK>
        end
        clear BlockInd
    end
    
    t=xlabel('Cortical depth');
    set(t,'fontsize',Fontsize);
    
    if exist('P', 'var')
        Sig = [];
        if exist('BOLD_Thresholds', 'var')
            if P(1,CondInd)<BOLD_Thresholds(1,CondInd)
                Sig = [Sig 'L ; p = ' sprintf('%.3f',(P(1,CondInd))) '\n']; %#ok<*AGROW>
            end
            %         if P(2,CondInd)<0.05
            %             Sig = [Sig 'Q ; p = ' sprintf('%.3f',(P(2,CondInd))) '\n'];
            %         end
            if P(2,CondInd)<BOLD_Thresholds(2,CondInd)
                Sig = [Sig 'C ; p = ' sprintf('%.3f',(P(2,CondInd)))];
            end
        else
            if P(1,CondInd)<0.05
                Sig = [Sig 'L ; p = ' sprintf('%.3f',(P(1,CondInd))) '\n']; %#ok<*AGROW>
            end
            %         if P(2,CondInd)<0.05
            %             Sig = [Sig 'Q ; p = ' sprintf('%.3f',(P(2,CondInd))) '\n'];
            %         end
            if P(2,CondInd)<0.05
                Sig = [Sig 'C ; p = ' sprintf('%.3f',(P(2,CondInd)))];
            end
        end
        
        if ~isempty(Sig)
            
            if NbConditions==4
                text(NbLayers/2,MAX-.5,sprintf(Sig))
            elseif any(CondInd==[4 8:12])
                text(NbLayers/2,Max-.1, sprintf(Sig))
            else
                text(NbLayers/2,MAX-.5,sprintf(Sig))
            end
            
        end
        clear Sig
    end
    
    set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
        'xticklabel', round(linspace(0,100,NbLayers)), 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    
    axis([1 NbLayers MIN MAX])
end



%%
if NbConditions==2
    subplot(m,n,1)
    %     t=ylabel('Main effect of attention on differentials');
    t=title('AV-A');
    set(t,'fontsize',Fontsize);
    axis([0.5 NbLayers+.5 MIN MAX])
    
    subplot(m,n,2)
    %     t=ylabel('Main effect of stimulation on differentials');
    t=title('AV-V');
    set(t,'fontsize',Fontsize+2);
    axis([0.5 NbLayers+.5 MIN MAX])
    
else
    subplot(m,n,1)
    t=ylabel('A attention');
    set(t,'fontsize',Fontsize+2);
    t=title('A stimulation');
    set(t,'fontsize',Fontsize);
    
    subplot(m,n,2)
    t=title('V stimulation');
    set(t,'fontsize',Fontsize);
    
    
    
    if NbConditions==6
        
        subplot(m,n,3)
        t=title('AV stimulation');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,4)
        t=ylabel('V attention');
        set(t,'fontsize',Fontsize+2);
        
    elseif NbConditions==4
        subplot(m,n,1)
        t=title('AV-A');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,2)
        t=title('AV-V');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,3)
        t=ylabel('V attention');
        set(t,'fontsize',Fontsize+2);
        
    elseif NbConditions==10
        subplot(m,n,3)
        t=title('AV stimulation');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,4)
        t=title('AV-A');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,5)
        t=title('AV-V');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,6)
        t=ylabel('V attention');
        set(t,'fontsize',Fontsize+2);
    else
        subplot(m,n,3)
        t=title('AV stimulation');
        set(t,'fontsize',Fontsize);
        
        subplot(m,n,4)
        t=title('AV - (A+V)');
        set(t,'fontsize',10);
        axis([0.5 NbLayers+.5 Min Max])
        
        subplot(m,n,5)
        t=ylabel('V attention');
        set(t,'fontsize',Fontsize+2);
        
        subplot(m,n,9)
        if strcmp(Name(1),'V');
            t=ylabel('V attention - A attention');
        else
            t=ylabel('A attention - V attention');
        end
        set(t,'fontsize',Fontsize+2);
        
        subplot(m,n,8)
        axis([0.5 NbLayers+.5 Min Max])
        
        for Ind=9:12
            subplot(m,n,Ind)
            t=xlabel('Cortical depth');
            set(t,'fontsize',Fontsize);
            axis([0.5 NbLayers+.5 Min Max])
        end
        clear Ind
        
        subplot(m,n,12)
        t=text(1.2,Max-1,['Nb Subject = ' num2str(length(DATA.Include))]);
        set(t,'fontsize',Fontsize);
        
        
        if isfield(DATA, 'VoxelCount')
            for iLayer=1:size(DATA.VoxelCount,2)
                t=text(iLayer-.25,Max-3,num2str(round(mean(DATA.VoxelCount(:,iLayer))/100)));
                set(t,'fontsize',Fontsize-2);
                t=text(iLayer-.25,Max-5,num2str(round(std(DATA.VoxelCount(:,iLayer))/100)));
                set(t,'fontsize',Fontsize-2);
            end
        end
        
    end
    
    if MainEffect
        subplot(m,n,9)
        t=ylabel('Main effect');
        set(t,'fontsize',Fontsize+2);
        
        subplot(m,n,4)
        t=title('Main effect');
        set(t,'fontsize',10);
    end
    
end

mtit(strrep(Name, '_', '-'), 'xoff', 0, 'yoff', 0.015);

%%
if Legend.Bool
    subplot(m,n,NbConditions)
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
            plot2svg(strcat(Name,'.svg'),gcf)
        case 3
            print(gcf, strcat(Name,'.pdf'), '-dpdf')
        case 4
            print(gcf, strcat(Name,'.tif'), '-dtiff')
            print(gcf, strcat(Name,'.pdf'), '-dpdf')
    end
end

end
