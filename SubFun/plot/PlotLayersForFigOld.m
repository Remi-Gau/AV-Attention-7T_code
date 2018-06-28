function PlotLayersForFig(DATA, Name, Visible)

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

Fontsize = 12;
Transparent = 0;


% Color for Conditions
COLOR =   [...
    255 150 150; ...
    150 255 150; ...
    150 220 255; ...
    255 0 0; ...
    0 255 0; ...
    0 0 255; ...
    150 150 150
    0 0 0];
COLOR = repmat([0 0 0],8,1);
COLOR=COLOR/255;

YLabel = {...
    'A attention';...
    '';...
    '';...
    'V attention';...
    '';...
    '';...
    '';...
    ''};

Title = {...
    'A stimulation';...
    'V stimulation';...
    'AV stimulation';...
    '';...
    '';...
    '';...
    'Main effect of MSI: $\frac{[AV-(A+V)]_{A}+[AV-(A+V)]_{V}}{2}$';...$\frac{[AV-(A+V)]_{A}+[AV-(A+V)]_{V}}{2}$ ([AV-(A+V)]_A + [AV-(A+V)]_V)/2
    'Main effect of attention: $\frac{[Att_V-Att_A]_{V}+[Att_V-Att_A]_{A}+[Att_V-Att_A]_{AV}}{3}$'}; %$\frac{[Att_V-Att_A]_{V}+[Att_V-Att_A]_{A}+[Att_V-Att_A]_{AV}}{3}$' ([Att_V-Att_A]_V + [Att_V-Att_A]_A + [Att_V-Att_A]_{AV})/3'


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

% COLOR_Subject = repmat([128 128 128],13,1);
COLOR_Subject=COLOR_Subject/255;


%%
Name = strrep(Name, ' ', '-');

if isfield(DATA, 'shadedErrorBar')
    Mean = DATA.shadedErrorBar.mean;
    ErrorBar = DATA.shadedErrorBar.errorbar;
    NbLayers = size(Mean,1);
    NbConditions = size(Mean,2);
end

if isfield(DATA, 'Subjects')
    Subjects = DATA.Subjects.Data;
    NbSubjects = size(Subjects,3);
    Transparent = 1;
    COLOR_Subject = COLOR_Subject(DATA.Include,:);
end

Beta = DATA.Beta;
if isfield(DATA, 'Thresholds')
    Thresholds = DATA.Thresholds;
end
if isfield(DATA, 'AllNullDist')
    AllNullDist = DATA.AllNullDist;
end

WithQuad = DATA.WithQuad;


m=4; n=5;
SubPlot = {1,2,3,6,7,8,[4 5 9 10], [11 12 16 17]};

Scatter = linspace(0,.4,11);


%%
figure('Name', Name, 'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);

box off

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

for CondInd=1:NbConditions
    
    if CondInd<7
        MAX = DATA.Max(1);
        MIN = DATA.Min(1);
    elseif CondInd==7
        MAX = DATA.Max(2);
        MIN = DATA.Min(2);
    elseif CondInd==8
        MAX = DATA.Max(3);
        MIN = DATA.Min(3);
    end
    
    
    subplot(m,n,SubPlot{CondInd})
    
    hold on; grid on;
    
    if CondInd<7
        shadedErrorBar(1:NbLayers, Mean(:,CondInd),ErrorBar(:,CondInd), ...
            {'LineWidth', 2, 'Color', COLOR(CondInd,:)}, Transparent)
        for SubjInd = 1:NbSubjects
            plot(1:NbLayers, flipud(Subjects(:,CondInd,SubjInd)), '--', ...
                'LineWidth', .5, 'Color', COLOR_Subject(SubjInd,:));
        end
    else
        shadedErrorBar(1:NbLayers, Mean(:,CondInd),ErrorBar(:,CondInd), ...
            {'Marker', '.', 'MarkerSize', 15, 'LineWidth', 3, 'Color', COLOR(CondInd,:)}, Transparent)
        for SubjInd = 1:NbSubjects
            plot(1:NbLayers, flipud(Subjects(:,CondInd,SubjInd)), '--', ...
                'LineWidth', 2, 'Color', COLOR_Subject(SubjInd,:));
            
            x = (1:.1:NbLayers)-mean(1:.1:NbLayers);
            b=squeeze(Beta(:,CondInd,SubjInd));
            if WithQuad
                y=b(1)*x+b(2)*x.^2+b(3);
            else
                y=b(1)*x+b(2);
            end
            plot(x-min(x)+1,y,':','LineWidth', 2, 'Color', COLOR_Subject(SubjInd,:))
            clear b
            
        end
        
        x = (1:.1:NbLayers)-mean(1:.1:NbLayers);
        b=mean(squeeze(Beta(:,CondInd,:))');
        if WithQuad
            y=b(1)*x+b(2)*x.^2+b(3);
        else
            y=b(1)*x+b(2);
        end
        plot(x-min(x)+1,y,':k','LineWidth', 3)
        clear b
        
    end
    
    plot([1 NbLayers], [0 0], '-k', 'LineWidth', 1)
    
    t=ylabel(YLabel{CondInd});
    set(t,'fontsize',Fontsize);
    
    if CondInd>6
        t=title(Title{CondInd}, 'Interpreter', 'Latex');
        set(t,'fontsize',Fontsize+5);
    else
        t=title(Title{CondInd});
        set(t,'fontsize',Fontsize);
    end
    
    
    if CondInd>6
        set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
            'xticklabel', linspace(0,100,NbLayers), 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
        t=xlabel('Cortical depth');
        set(t,'fontsize',Fontsize);
    else
        set(gca,'tickdir', 'out', 'xtick', 1:NbLayers , ...
            'xticklabel', [], 'ticklength', [0.01 0.01], 'fontsize', Fontsize)
    end
    
    axis([1 NbLayers MIN MAX])
    
    
    % Inset with betas
    if CondInd>=7
        
        tmp = squeeze(Beta(:,CondInd,:));
        
        for i=1:size(tmp,1)
            
            Lim = ceil(max(abs([min(tmp(i,:))-.02 max(tmp(i,:))+.02]))*100)/100;
            
            if CondInd==7
                XY = [0.65 0.35];
            else
                XY = [0.465 0.14];
            end
            
            axes('Position',[XY(1)+0.1*(i-1) XY(2) .04 .09])
            box off; hold on;
            
            if exist('Thresholds','var')
                [H(i),P(i)] = ttest(tmp(i,:), 0, 'alpha', Thresholds(i,CondInd-6));
            end
            
            if exist('AllNullDist','var')
                NullDist = AllNullDist(:,i,CondInd-6);
                
                P(i) = sum(NullDist>abs(mean(tmp(i,:))))/numel(NullDist);
                if P(i)<.05
                    H(i) = 1;
                else
                    H(i) = 0;
                end
                
                distributionPlot(NullDist, 'xValues', 0.5, 'histOri', 'left', ...
                    'histOpt', 0, 'divFactor', 100, 'showMM', 0, 'distWidth', .5 );
                
                plot([0.5 0.9],[NullDist(floor(0.95*numel(NullDist))) NullDist(floor(0.95*numel(NullDist)))],'-r')
                
            end
            
            for SubjInd=1:size(tmp,2)
                plot(1.2+Scatter(SubjInd), tmp(i,SubjInd), ...
                    'Marker', '.', 'MarkerEdgeColor', COLOR_Subject(SubjInd,:), ...
                    'MarkerFaceColor', COLOR_Subject(SubjInd,:), 'MarkerSize', 15)
            end
            
            plot([0 1.8], [0 0], ':k', 'LineWidth', .5)
            
            if exist('AllNullDist','var')
                h = errorbar(1,abs(mean(tmp(i,:),2)),nansem(abs(tmp(i,:)),2), '.k');
            else
            h = errorbar(1,mean(tmp(i,:),2),nansem(tmp(i,:),2), '.k');
            end
            set(h, 'linewidth', 1, 'MarkerSize', 10)
            
            Sig = [];
            
            
            if P(i)<0.001
                Sig = sprintf('ES=%.3f \np<0.001 ',...
                    abs(mean(tmp(i,:))/std(tmp(i,:))));
            else
                Sig = sprintf('ES=%.3f \np=%.3f ',...
                    abs(mean(tmp(i,:))/std(tmp(i,:))), P(i));
            end
            
            t = text(1,Lim+Lim*50/100,sprintf(Sig));
            set(t,'fontsize',Fontsize-2);
            
            if H(i)==1
                set(t,'color','r');
            end
            
            clear Sig
            
            switch i
                case 1
                    xTickLabel = 'L';
                case 2
                    if size(tmp,1)==2
                        xTickLabel = 'C';
                    else
                        xTickLabel = 'Q';
                    end
                case 3
                    xTickLabel = 'C';
            end
            
            
            set(gca,'tickdir', 'out', 'xtick', 1.1:3.1 , 'xticklabel', xTickLabel, ...
                'ytick', linspace(Lim*-1, Lim, 5) , 'yticklabel', linspace(Lim*-1, Lim, 5), 'ticklength', [0.05 0.05], 'fontsize', Fontsize-2)
            t=ylabel('betas');
            set(t,'fontsize',Fontsize-3);
            
            axis([0 1.8 Lim*-1 Lim])
            
        end
        
    end
    
end

mtit(strrep(Name(1:end-4),'_',' '), 'xoff', 0, 'yoff', +0.03, ...
    'fontsize', 16)


% plot2svg(strcat(Name, '_', num2str(NbLayers), 'Layers.svg'),gcf)
% print(gcf, strcat(Name,'_', num2str(NbLayers), 'Layers.tif'), '-dtiff')
print(gcf, strcat(Name, '_',num2str(NbLayers), 'Layers.pdf'), '-dpdf')


end
