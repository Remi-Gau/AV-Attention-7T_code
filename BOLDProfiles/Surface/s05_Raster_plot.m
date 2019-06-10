%%
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

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

ToPlot={'Cst','Lin'};

% ROI_Groups = {1:4;5:8;9:12;13:16;17:20;21:24;25:28};
% ROI_Groups_names = {'A1-PT-V123';'V123-ActDeact';'A1-PT-A-ActDeact';'V123-A-ActDeact';...
%     'A1-PT-V-ActDeact';'V123-V-ActDeact';'A1-PT-ActDeact'};
ROI_Groups = {1:4};
ROI_Groups_names = {'A1-PT-V123'};

load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
% ToPermute = [];


%% Grp level  ; raster stim = f(Percentile of V)
close all
Cdt =[2 1;2 2;2 3];
Name={'A','V','AV'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                tic
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                parfor iSubj = 1:size(All_Profiles_V,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-4 4])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI==ROI_Groups{iRoiGrp}(1)
                    title([Name{Cdt(iCdt,2)} ' stimulus'])
                end
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .15, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s %s stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, [-4 4])
                end
                
                clear X
                
                toc
                
            end
            
        end
        
        mtit(['Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_V_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
    
end


%% Grp level  ; raster stim = f(Percentile of A)
Cdt =[1 1;1 2;1 3];
Name={'A','V','AV'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:3
            
            for iROI = ROI_Groups{iRoiGrp}
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                clear X_sort Profiles Sorting_Raster
                
                parfor iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),3,iCdt+3*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-4 4])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .15, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s %s stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, [-4 4])
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Percentile A stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Baseline', ...
            ['GrpLvl_Raster_Baseline_A_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster Att = f(Percentile of V)
close all
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:4
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles Sorting_Raster
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_AttfV{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                
                subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch{iCdt});
                    set(t,'fontsize',6)
                else
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .1, .1, .04, .04, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf(sprintf('%s\nPerc %s V stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}));
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, [-2 2])
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile V stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
%         print(gcf, fullfile(FigureFolder,'Attention', ...
%             ['GrpLvl_Raster_Attention_V_stim_' ToPlot{iToPlot} '_' ...
%             ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
    end
end


%% Grp level  ; raster Att = f(Percentile of A)
Name={'Att_A - Att_V', 'A stim ; Att_A - Att_V', 'V stim ; Att_A - Att_V', 'AV stim ; Att_A - Att_V'};
NameSwitch={'Att_V - Att_A', 'A stim ; Att_V - Att_A', 'V stim ; Att_V - Att_A', 'AV stim ; Att_V - Att_A'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:4
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles Sorting_Raster
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_Profiles_A,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_AttfA{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),4,iCdt+4*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch{iCdt});
                    set(t,'fontsize',6)
                else
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .1, .1, .04, .04, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s A stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, [-2 2])
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile A stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'Attention', ...
            ['GrpLvl_Raster_Attention_A_stim_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
        
    end
end



%% Grp level  ; raster CrossMod = f(Percentile of V stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,1,iROI};
                    X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI == ROI_Groups{iRoiGrp}(1)
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                YLabel = sprintf('%s\nPerc %s V stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                
                clear X
                
            end
            
        end
        
        mtit(['CrossMod = f(Percentile V Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_CrossMod_Stim_V_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
    end
end


%% Grp level  ; raster CrossMod = f(Percentile of A stim)
close all
Name = {'AV-A';'AV-V'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-2 2])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if iROI == ROI_Groups{iRoiGrp}(1)
                    t=title(Name{iCdt});
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                YLabel = sprintf('%s\nPerc %s A stim',...
                    strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot});
                PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                
                clear X
                
            end
            
        end
        
        mtit(['CrossMod = f(Percentile A Stim) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_CrossMod_Stim_A_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
        
    end
end


%% Grp level  ; raster Att = f(Percentile of CrossMod)
close all
Cdt =[1 4;2 4];
Name = {'AV-A';'AV-V'};
NameAtt={'Att_A - Att_V'};
NameSwitch={'Att_V - Att_A'};
ToPlot={'Cst','Lin'};

for iToPlot = 1:numel(ToPlot)
    
    for iRoiGrp = 1:numel(ROI_Groups)
        
        figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
        
        for iCdt = 1:2
            
            for iROI = ROI_Groups{iRoiGrp}
                
                clear X_sort Profiles
                
                NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
                
                for iSubj = 1:size(All_X_sort_Cross,1)
                    X_sort(iSubj,:) = All_X_sort_Att_Cross{iSubj,iToPlot,iCdt,iROI};
                    Profiles(:,:,iSubj) = All_Profiles_Att_Cross{iSubj,iToPlot,iCdt,iROI};
                end
                
                X_sort = X_sort(Subj2Include,:);
                Profiles = Profiles(:,:,Subj2Include);
                
                [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
                
                subplot(numel(ROI_Groups{iRoiGrp}),2,iCdt+2*(iROI-ROI_Groups{iRoiGrp}(1)))
                hold on
                
                colormap(ColorMap);
                imagesc(imgaussfilt(mean(Profiles,3),[20 .001]), [-1 1])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 4)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',6)
                
                if ROI(iROI).name(1)=='V'
                    t=title(NameSwitch);
                    set(t,'fontsize',6)
                else
                    t=title(NameAtt);
                    set(t,'fontsize',6)
                end
                
                
                ax = gca;
                
                rho = atanh(rho);
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, .275, .1, .05, .05, [0.9 1.4 -1 1], ToPermute);
                
                if iCdt==1
                    YLabel = sprintf('%s\nPerc %s',...
                        strrep(ROI(iROI).name,'_','-'), Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                else
                    YLabel = sprintf('Perc %s',...
                        Name{Cdt(iCdt,1)});
                    PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1)
                    
                end
                
                clear X
                
            end
            
        end
        
        mtit(['Att = f(Percentile CrossMod) - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
        
        print(gcf, fullfile(FigureFolder,'CrossSensory', ...
            ['GrpLvl_Raster_Att_CrossMod_' ToPlot{iToPlot} '_' ...
            ROI_Groups_names{iRoiGrp} '_CV.tif']), '-dtiff')
        
        
    end
    
end

