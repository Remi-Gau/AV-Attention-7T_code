clear
clc

% Color map
X = 0:0.001:1;
R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);
ColorMap = [R' G' B'];
clear R G B

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

StartFolder='/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

load(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster_CV.mat'))

Subj2Include = true(11,1);
% Subj2Include([4 11]) = false;

for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
% ToPermute = [];

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

Color(10,:) = [1,0,0];

ToPlot={'Cst','Lin'};

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

CorticalDepth = round(linspace(100,0,NbLayers+2));
CorticalDepth([1 end]) = [];

% ROI_Groups = {1:4;5:8;9:12;13:16;17:20;21:24;25:28};
% ROI_Groups_names = {'A1-PT-V123';'V123-ActDeact';'A1-PT-A-ActDeact';'V123-A-ActDeact';...
%     'A1-PT-V-ActDeact';'V123-V-ActDeact';'A1-PT-ActDeact'};

ROI_Groups = {1:4};
ROI_Groups_names = {'A1-PT-V123'};


load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')



%% Grp level  ; raster AV-A = f(Percentile of A stim)
close all
SubPlot = [1 3 2 4];
CLIM = [-3 3];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI}; %#ok<*SAGROW>
                    Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(2,2,SubPlot(iROI))
                hold on
                
                colormap(ColorMap);
                imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
              
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s A stim',ROI(iROI).name,ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .25, .25, .08, .08, [0.9 1.4 -2 1], ToPermute)
                
                if iROI>2
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['AV-A = f(Percentile A Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_AV-A_Stim_A_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
    
end


%% Grp level  ; raster AV-V = f(Percentile of V stim)
close all
SubPlot = [1 3 2 4];
CLIM = [-3 3];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                    X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(2,2,SubPlot(iROI))
                hold on
                
                colormap(ColorMap);
                imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                ax = gca;
                
                YLabel = sprintf('%s\nPerc %s V stim',ROI(iROI).name,ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .25, .25, .08, .08, [0.9 1.4 -1.1 0.5], ToPermute)
                
                if iROI>2
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['AV-V = f(Percentile V Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_AV-V_Stim_V_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
        
    end
    
end


%% Grp level  ; raster A stim = f(Percentile of V)
close all
Cdt =[2 1];
Name={'A','V'};
SubPlot = [1 3 2 4];
CLIM = [-3 3];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_Profiles_V,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                    X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(2,2,SubPlot(iROI))
                hold on
                
                colormap(ColorMap);
                %             imagesc(mean(fliplr(Profiles),3), [-3 3])
                imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                title(ROI(iROI).name)
                
                ax = gca;
                
                if iROI<3
                    YLabel = sprintf('Perc %s %s stim', ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                else
                    YLabel = [];
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 0, Sorting_Raster, CLIM)
                end
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .23, .23, .08, .08, [0.9 1.4 -0.25 0.75], ToPermute)
                
                if iROI>2
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['A stim = f(Percentile V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.03)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_V_stim_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster V stim = f(Percentile of A)
close all
Cdt =[1 2;1 2];
Name={'A','V'};
SubPlot = [1 3 2 4];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                
                subplot(2,2,SubPlot(iROI))
                hold on
                
                colormap(ColorMap);
                if iROI<3
                    imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-.8 .8])
                else
                    imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-4 4])
                end
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                title(ROI(iROI).name)
                
                ax = gca;
                
                if iROI<3
                    YLabel = sprintf('Perc %s %s stim', ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                else
                    YLabel = [];
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 0, Sorting_Raster, CLIM)
                end
                
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .23, .23, .08, .08, [0.9 1.4 -.75 .75], ToPermute)
                
                if iROI<3
                    PlotColorBar(ax, ColorMap, [-.8 .8])
                else
                    PlotColorBar(ax, ColorMap, [-4 4])
                end
                
                clear X
                
            end
            
        end
        
        mtit(['V stim = f(Percentile A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_A_stim_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster V/A stim = f(Percentile of A/V)
Cdt =[1 2;2 1];
Name={'A','V'};
SubPlot = [1 3 2 4];

CLIM = [-1.2 1.2];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        
        for iROI = ROI_Groups{iRoiGrp}
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            if iROI<3
                iCdt = 2;
            else
                iCdt = 1;
            end
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_A,1)
                if iROI<3
                    Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
                else
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                    X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
                end
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            subplot(2,2,SubPlot(iROI))
            hold on
            
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            if iROI<3
                t=title(sprintf('%s\nV stim = f(Percentile A stim)',ROI(iROI).name));
            else
                t=title(sprintf('%s\nA stim = f(Percentile V stim)',ROI(iROI).name));
            end
            
            ax = gca;
            
            YLabel = sprintf('Perc %s %s stim', ToPlot{iToPlot}, Name{Cdt(iCdt,2)});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
            
            rho = atanh(rho);
            PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .23, .23, .08, .08, [0.9 1.4 -.5 .5], ToPermute)
            
            if iROI>2
                PlotColorBar(ax, ColorMap, CLIM)
            end
            
            
            clear X
            
        end
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_A_V_stim_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster Att = f(Percentile of V)
close all
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

CLIM = [-2 2];

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                    X_sort(iSubj,:) = All_X_sort_AttfV{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-1))
                hold on
                
                colormap(ColorMap);
                imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                if iROI==1
                    title(Name{iCdt})
                    
                elseif iROI==3
                    title(NameSwitch{iCdt})
                    
                elseif iROI==numel(ROI)
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', [],'yticklabel', [], ...
                        'ticklength', [0.01 0], 'fontsize', 10)
                    t=xlabel('cortical depth');
                    set(t,'fontsize',10)
                end
                
                ax = gca;
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s V stim', ROI(iROI).name,ToPlot{iToPlot});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                end
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .135, .075, .06, .06, [0.9 1.4 -.5 .5], ToPermute)
                
                if iCdt==3
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['GrpLvl_Raster_Attention_V_stim_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster Att = f(Percentile of A)
close all
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_AttfA{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-1))
                hold on
                
                colormap(ColorMap);
                imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-1.5 1.5])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                if iROI==1
                    title(Name{iCdt})
                    
                elseif iROI==3
                    title(NameSwitch{iCdt})
                    
                elseif iROI==numel(ROI)
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', [],'yticklabel', [], ...
                        'ticklength', [0.01 0], 'fontsize', 10)
                    t=xlabel('cortical depth');
                    set(t,'fontsize',10)
                end
                
                ax = gca;
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s A stim', ROI(iROI).name,ToPlot{iToPlot});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)                    
                end
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .135, .075, .06, .06, [0.9 1.4 -.5 .5], ToPermute)
                
                if iCdt==3
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['GrpLvl_Raster_Attention_A_stim_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
        
    end
end


%% Grp level  ; raster CrossMod = f(Percentile of V stim)
close all
Name = {'AV-A';'AV-V'};

CLIM = [-0.5 0.5];

iRoiGrp = 1;

for iToPlot = 2
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1
        
        for iROI = 1:2
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            subplot(2,1,iROI)
            hold on
            
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel', CorticalDepth, ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==numel(ROI)
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            ax = gca;
            
            YLabel = sprintf('%s\nPerc %s V stim',ROI(iROI).name,ToPlot{iToPlot});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
            
            PlotColorBar(ax, ColorMap, CLIM)
            
            rho = atanh(rho);
            PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .68, .25, .08, .08, [0.9 1.4 -.081 .081], ToPermute)
            
            
            
            
            clear X
            
        end
        
    end
    
    mtit(['CrossMod = f(Percentile V Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Raster_CrossMod_Stim_V_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
    
    
end


%% Grp level  ; raster CrossMod = f(Percentile of A stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 2
        
        for iROI = 3:4
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            subplot(2,1,iROI-2)
            hold on
            
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-.5 .5])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==numel(ROI)
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            ax = gca;
            
            YLabel = sprintf('%s\nPerc %s A stim',ROI(iROI).name,ToPlot{iToPlot});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
            
            rho = atanh(rho);
            PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .68, .25, .08, .08, [0.9 1.4 -.5 .5], ToPermute)
            
            PlotColorBar(ax, ColorMap, [-1 1])
            
            clear X
            
        end
        
    end
    
    mtit(['CrossMod = f(Percentile A Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Raster_CrossMod_Stim_A_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
    
    
end


%% Grp level  ; raster Att = f(Percentile of CrossMod)
close all
Cdt =[1 4;2 4];
Name = {'AV-A';'AV-V'};
NameAtt={'Att_A - Att_V'};
NameSwitch={'Att_V - Att_A'};
ToPlot={'Cst','Lin'};
SubPlot = [1 3 2 4];

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    iSubplot = 1;
    
    for iCdt = 1:2
        
        switch iCdt
            case 1
                RoiToPlot = [1 2];
            case 2
                RoiToPlot = [3 4];
        end
        
        for iROI = RoiToPlot
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_Cross{iSubj,iToPlot,iCdt,iROI};
                X_sort(iSubj,:) = All_X_sort_Att_Cross{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Att_Cross{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            subplot(2,2,SubPlot(iSubplot))
            iSubplot = iSubplot + 1;
            hold on
            
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), [-0.8 .8])
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            if iROI==1
                title(NameAtt{1})
                
            elseif iROI==3
                title(NameSwitch{1})
            end
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            ax = gca;
            
            YLabel = sprintf('%s\nPerc %s',ROI(iROI).name, Name{Cdt(iCdt,1)});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
            
            rho = atanh(rho);
            PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .25, .255, .08, .08, [0.9 1.4 -.5 .5], ToPermute)
            
            if iROI>2
                PlotColorBar(ax, ColorMap, [-0.8 .8])
            end
            
            clear X
            
        end
        
    end
    
    mtit(['Att = f(Percentile CrossMod) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Raster_Att_CrossMod_Final_' ToPlot{iToPlot} '_' ...
                    ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
    
    
end




return

%% Grp level  ; raster CrossMod = f(Percentile of V stim)
close all
Name = {'AV-A';'AV-V'};

CLIM = [-0.5 0.5];

iRoiGrp = 1;

for iToPlot = 2
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1
        
        for iROI = 1
            
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            subplot(2,1,iROI)
            hold on
            
            colormap(ColorMap);
            imagesc(mean(imgaussfilt(Profiles,[20 .001]),3), CLIM)
            
            axis([0.5 6.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel', CorticalDepth, ...
                'ytick', [],'yticklabel', [], ...
                'ticklength', [0.01 0], 'fontsize', 10)
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==numel(ROI)
                set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            ax = gca;
            
            YLabel = sprintf('%s\nPerc %s V stim',ROI(iROI).name,ToPlot{iToPlot});
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
            
            PlotColorBar(ax, ColorMap, CLIM)
            
            
            subplot(2,1,2)
            ax = gca;
            axis off
            rho = atanh(rho);
            PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0.1, 0, .5, .4, [0.9 1.4 -.081 .081], ToPermute)
            
            
            
            
            clear X
            
        end
        
    end
    
    mtit(['CrossMod = f(Percentile V Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Raster_CrossMod_Stim_V_Final_' ToPlot{iToPlot} '_' ...
                    'A1_CV.tif']), '-dtiff')
    
    
end




return

%% Grp level  ; profiles A stim = f(Percentile of V)
SubPlot = [1 3 2 4];

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1
        
        for iROI = 1:numel(ROI)
            
            clear X_sort Profiles
            
            for iSubj = 1:size(All_Profiles_V,1)
                Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
                X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            X_sort = X_sort(Subj2Include,:);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
                X_sort_Quant(:,iPerc-1,:) = mean(X_sort(:,IdxToAvg((iPerc-1):iPerc),:),2);
                
                Y = squeeze(Profiles_Quant(iPerc-1,:,:));
                
                X = [];
                for iSubj=1:size(Y,2)
                    X((1:6)+6*(iSubj-1),(1:2)+2*(iSubj-1)) = DesMat; 
                end
                
                Y = Y(:);
                B = pinv(X)*Y;
                
                Cst_tmp = B(1:2:size(X,2),:);
                Lin_tmp = B(2:2:size(X,2),:);
                
                if iToPlot==1
                    Y_sort_Quant(:,iPerc-1,:)=Cst_tmp;
                else
                    Y_sort_Quant(:,iPerc-1,:)=Lin_tmp;
                end
                
            end
            
            for iSubj=1:size(Y_sort_Quant,1)
                R=corrcoef(X_sort_Quant(iSubj,:),Y_sort_Quant(iSubj,:));
                rho(iSubj) = R(1,2);
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(2,2,SubPlot(iROI))
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            if iROI<3
                axis([0.5 6.5 -.5 4.5])
            else
                axis([0.5 6.5 -4.5 0.5])
            end
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                'ticklength', [0.01 0], 'fontsize', 8)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            title(ROI(iROI).name)
            
            t=ylabel('Param est.');
            set(t,'fontsize',10)
            
        end
        
    end
    
    mtit(['A stim profiles = f(top/bottom deciles V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'Baseline', ...
        ['GrpLvl_Profiles_Baseline_V_stim_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end


%% Grp level  ; profiles V stim = f(Percentile of A)
SubPlot = [1 3 2 4];

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 2
        
        for iROI = 1:numel(ROI)
            
            clear X_sort Profiles Profiles_Quant
            
            for iSubj = 1:size(All_Profiles_A,1)
                Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
                X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            X_sort = X_sort(Subj2Include,:);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
                X_sort_Quant(:,iPerc-1,:) = mean(X_sort(:,IdxToAvg((iPerc-1):iPerc),:),2);
                
                Y = squeeze(Profiles_Quant(iPerc-1,:,:));
                
                X = [];
                for iSubj=1:size(Y,2)
                    X((1:6)+6*(iSubj-1),(1:2)+2*(iSubj-1)) = DesMat; 
                end
                
                Y = Y(:);
                B = pinv(X)*Y;
                
                Cst_tmp = B(1:2:size(X,2),:);
                Lin_tmp = B(2:2:size(X,2),:);
                
                if iToPlot==1
                    Y_sort_Quant(:,iPerc-1,:)=Cst_tmp;
                else
                    Y_sort_Quant(:,iPerc-1,:)=Lin_tmp;
                end
                
            end
            
            for iSubj=1:size(Y_sort_Quant,1)
                R=corrcoef(X_sort_Quant(iSubj,:),Y_sort_Quant(iSubj,:));
                rho(iSubj) = R(1,2);
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(2,2,SubPlot(iROI))
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            if iROI<3
                axis([0.5 6.5 -1.5 1])
            else
                axis([0.5 6.5 -.5 7.5])
            end
            
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                'ticklength', [0.01 0], 'fontsize', 8)
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            title(ROI(iROI).name)
            
            t=ylabel('Param est.');
            set(t,'fontsize',10)
            
        end
        
        
        
    end
    
    mtit(['V stim profiles = f(top/bottom deciles A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'Baseline', ...
        ['GrpLvl_Profiles_Baseline_A_stim_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end


%% Grp level  ; Profiles Att = f(Percentile of V)
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1:3
        
        for iROI = 1:numel(ROI)
            
            clear X_sort Profiles
            
            for iSubj = 1:size(All_Profiles_AttfV,1)
                Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(numel(ROI),3,iCdt+3*(iROI-1))
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            
            switch iCdt
                case 1
                    axis([0.5 6.5 -2 3.5])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
                    
                case 2
                    axis([0.5 6.5 -2 3])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
                case 3
                    axis([0.5 6.5 -3 4.5])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
            end
            
            
            
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==3
                title(NameSwitch{iCdt})
                
            elseif iROI==numel(ROI)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            
            if iCdt==1
                t=ylabel(sprintf('%s\nParam est.',...
                    ROI(iROI).name));
                set(t,'fontsize',10)
            end
            
        end
        
    end
    
    mtit(['profiles = (top/bottom deciles V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'Attention', ...
        ['GrpLvl_Profiles_Attention_V_stim_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end


%% Grp level  ; Profiles Att = f(Percentile of A)
close all
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1:3
        
        for iROI = 1:numel(ROI)
            
            clear X_sort Profiles
            
            for iSubj = 1:size(All_Profiles_AttfV,1)
                Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(numel(ROI),3,iCdt+3*(iROI-1))
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            
            switch iCdt
                case 1
                    axis([0.5 6.5 -1.5 2])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:.5:20,'yticklabel', -10:.5:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
                    
                case 2
                    axis([0.5 6.5 -2 2.5])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:1:20,'yticklabel', -10:1:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
                case 3
                    axis([0.5 6.5 -2 2])
                    set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                        'ytick', -10:.5:20,'yticklabel', -10:.5:20, ...
                        'ticklength', [0.01 0], 'fontsize', 8)
            end
            
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==3
                title(NameSwitch{iCdt})
                
            elseif iROI==numel(ROI)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            
            if iCdt==1
                t=ylabel(sprintf('%s\nParam est.',...
                    ROI(iROI).name));
                set(t,'fontsize',10)
            end
            
        end
        
    end
    
    mtit(['profiles = (top/bottom deciles A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'Attention', ...
        ['GrpLvl_Profiles_Attention_A_stim_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end


%% Grp level  ; Profiles CrossMod = f(Percentile of V Stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 1
        
        for iROI = 1:2 %numel(ROI)
            
            clear X_sort Profiles
            
            for iSubj = 1:size(All_Profiles_Cross,1)
                Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(2,1,iROI)
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            axis([0.5 6.5 -2 4.5])
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', -12:1:24,'yticklabel', -12:1:24, ...
                'ticklength', [0.01 0], 'fontsize', 8)
            
            
            if iROI==1
                title(Name{iCdt})
                
            elseif iROI==numel(ROI)
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
            end
            
            
            if iCdt==1
                t=ylabel(sprintf('%s\nParam est.',...
                    ROI(iROI).name));
                set(t,'fontsize',10)
            end
            
        end
        
    end
    
    mtit(['CrossMod = (top/bottom deciles Stim V) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Profiles_CrossMod_Stim_V_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end


%% Grp level  ; Profiles CrossMod = f(Percentile of A Stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
    
    for iCdt = 2
        
        for iROI = 3:4
            
            clear X_sort Profiles
            
            for iSubj = 1:size(All_Profiles_Cross,1)
                Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            Profiles = Profiles(:,:,Subj2Include);
            
            IdxToAvg = floor(linspace(1,size(Profiles,1),10+1));
            
            for iPerc = 2:numel(IdxToAvg)
                Profiles_Quant(iPerc-1,:,:) = mean(Profiles(IdxToAvg((iPerc-1):iPerc),:,:));
            end
            
            Profiles_Quant = fliplr(Profiles_Quant);
            
            subplot(2,1,iROI-2)
            hold on
            grid on
            plot([0 7], [0 0], ':k', 'LineWidth', 2)
            for iQuant = [1 size(Profiles_Quant,1)]
                shadedErrorBar(1:6, nanmean(Profiles_Quant(iQuant,:,:),3),nansem(Profiles_Quant(iQuant,:,:),3), ...
                    {'Marker', '.', 'MarkerSize', 10, 'LineStyle', '-', 'LineWidth', 1, 'Color', Color(iQuant,:)}, 'on')
            end
            
            axis([0.5 6.5 -2.2 2.5])
            set(gca,'tickdir', 'out', 'xtick', 1:6,'xticklabel',  CorticalDepth, ...
                'ytick', -12:1:24,'yticklabel', -12:1:24, ...
                'ticklength', [0.01 0], 'fontsize', 8)
            
            title(Name{iCdt})
            
            t=xlabel('cortical depth');
            set(t,'fontsize',10)
            
            t=ylabel(sprintf('%s\nParam est.',...
                ROI(iROI).name));
            set(t,'fontsize',10)
            
        end
        
    end
    
    mtit(['CrossMod = (top/bottom deciles Stim A) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.035)
    
    print(gcf, fullfile(FigureFolder,'CrossSensory', ...
        ['GrpLvl_Profiles_CrossMod_Stim_A_Final_' ToPlot{iToPlot} '.tiff']), '-dtiff')
    
end



