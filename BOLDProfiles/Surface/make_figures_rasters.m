%% make raster figures for the article

clear; close all; clc;

% CodeFolder = '/home/remi/github/AV-Attention-7T_code';
CodeFolder = 'D:\github\AV-Attention-7T_code';

% inputs (where the OSF data have been downloaded: https://osf.io/63dba/)
% DataFolder = '/home/remi/';
DataFolder = 'D:\';
DataFolder = fullfile(DataFolder, 'Dropbox', 'PhD', 'Experiments', 'AV_Integration_7T');

Results_Folder = fullfile(DataFolder, 'DataToExport');

% output folder
FigureFolder = fullfile(CodeFolder, 'Figures');
mkdir(FigureFolder)


%% data files parameters
NbLayers = 6;

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '15';...
    '16'
    ];

ToPlot={'Cst','Lin'};


%% Figures parameters

FigDim = [100, 100, 1000, 700];
Visibility = 'on';

%to pplot each subject's regresion coeficient, the group mean and the p value
plot_reg_coeff = 1;

% row demeaning of rasters plots for linear shape parameters
demean_raster = 1;


%% Design matrix for laminar GLM
DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)']; % in case we want a
% quadratic component
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


%% get things ready
NbSubj = size(SubjectList,1);
Subj2Include = true(NbSubj,1);

if plot_reg_coeff
    NbSubPlot = 4;
else
    NbSubPlot = 3;
end

% create permutations for exact sign permutation test
for iSubj=1:NbSubj
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:}); clear sets
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
clear a b c d e f g h i j k
ToPermute = []; %  comment to run permutation test

% add dependencies
addpath(genpath(fullfile(CodeFolder, 'SubFun')))
Get_dependencies('/home/remi')
Get_dependencies('D:')

% Color map
X = 0:0.001:1;
ColorMap = seismic(size(X,2));
clear X

% load data
load(fullfile(Results_Folder,'MinNbVert.mat'), 'MinVert')
load(fullfile(Results_Folder,'Raster_CV.mat'))



%% raster CrossMod = f(Percentile of V stim)
% for figure 3
% will create extra rasters not used in the end

close all

plot_reg_coeff = 0;
demean_raster = 0;

Sorting_Cdt = 2;
mkdir(fullfile(FigureFolder,'CrossSensory'));
Name = {'AV-A';'AV-V'};

if plot_reg_coeff
    NbSubPlot = 4;
else
    NbSubPlot = 3;
end

for iToPlot = 1:2
    
    for iCdt = 1:2
        
        for iROI = 1:4
            
            if iROI==1
                CLIM = [-0.15 0.15]; % for main raster
                CLIM2 = CLIM*1.1; % for sorting raster
            elseif iROI == 2
                CLIM = [-0.8 0.8];
                CLIM2 = CLIM;
            else
                CLIM = [-1 1];
                CLIM2 = [-1 1];
            end
            
            fig_name = ['GrpLvl_Raster_' Name{iCdt} ...
                '_fV_' ToPlot{iToPlot} '_CV_' ROI(iROI).name];
            
            disp(fig_name)
            
            figure('name', strrep(fig_name, '_', ' '), 'Position', ...
                FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            % get mininum number of bin for that raster (set by the subject
            % with the smallest number of vertices with data in his/her
            % ROI)
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:, :, iSubj) = All_Profiles_V{iSubj, iToPlot, Sorting_Cdt, iROI};
                X_sort(iSubj, :) = All_X_sort_Cross_V{iSubj, iToPlot, iCdt, iROI};
                Profiles(:, :, iSubj) = All_Profiles_Cross_V{iSubj, iToPlot, iCdt, iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            % computes vertex wise correlation / regression between predictor
            % and predicted rasters plots
            [rho,slope] = CorRegRaster(Profiles, DesMat, iToPlot, X_sort);
            
            
            %% plot main raster
            subplot(1, NbSubPlot, 2:3)
            PlotRectangle(NbLayers, 10, 1, 0)
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            subplot(1, NbSubPlot, 2:3)
            hold on
            colormap(ColorMap);
            
            % plot raster (add some gaussian smoothing across vertices)
            raster_to_plot = mean(imgaussfilt(Profiles,[size(Profiles,1)/100 .0001]), 3);
            if demean_raster && strcmp(ToPlot{iToPlot}, 'Lin')
                mnval = mean(raster_to_plot, 2);
                for i=1:size(raster_to_plot, 1)
                    raster_to_plot(i,:) = raster_to_plot(i,:) - mnval(i);
                end
            end
            imagesc(raster_to_plot, CLIM)
            axis([0.5 NbLayers+.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            
            %% plot sorting raster
            subplot(1, NbSubPlot, 1)
            PlotRectangle(NbLayers, 10, 1, 0)
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            subplot(1, NbSubPlot, 1)
            hold on
            colormap(ColorMap);
            
            % plot raster (add some gaussian smoothing across vertices)
            raster_to_plot = mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3);
            if demean_raster && strcmp(ToPlot{iToPlot}, 'Lin')
                mnval = mean(raster_to_plot, 2);
                for i=1:size(raster_to_plot, 1)
                    raster_to_plot(i,:) = raster_to_plot(i,:) - mnval(i);
                end
            end
            imagesc(raster_to_plot, CLIM2)
            axis([0.5 NbLayers+.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            
            %% plot sorting value
            ax = gca;
            YLabel = [];
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, 0)
            set(gca,'fontsize', 20)
            
            % plot regression coefficient
            if plot_reg_coeff
                subplot(1, NbSubPlot, 4) %#ok<*UNRCH>
                ax = gca;
                axis off
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
                    [0.9 1.6 -0.1 0.5], ToPermute)
                set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
            end
            
            
            %% saves table results and figure
            SavedTxt = fullfile(FigureFolder,'CrossSensory', [fig_name '.csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, [Name{iCdt} '_fV'], slope, ToPermute)
            
            
            print(gcf, fullfile(FigureFolder,'CrossSensory', [fig_name '.tif']), '-dtiff')
            
        end
        
    end
    
    
end


%% raster Cdt = f(Percentile of V stim)
% for extra figure to show that response to one condition does not predict
% responses to other conditions which would be expected if non-neurophysio
% (big vasculature) factors were most important.

close all

plot_reg_coeff = 0;
demean_raster = 0;

Sorting_Cdt = 2;

mkdir(fullfile(FigureFolder,'Baseline'));

Cdt = [1 2 3];
Name = {'A', 'V', 'AV'};

CLIM = [-3 3];
CLIM2 = [-0.8 .8];

if plot_reg_coeff
    NbSubPlot = 4;
else
    NbSubPlot = 3;
end

for iToPlot = 1:2
    
    for iCdt = 1:3
        
        for iROI = 1:4
            
            fig_name = ['GrpLvl_Raster_' Name{iCdt} ...
                '_fV_' ToPlot{iToPlot} '_CV_' ROI(iROI).name];
            
            disp(fig_name)
            
            figure('name', strrep(fig_name, '_', ' '), 'Position', ...
                FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            % get mininum number of bin for that raster (set by the subject
            % with the smallest number of vertices with data in his/her
            % ROI)
            NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
            
            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:, :, iSubj) = All_Profiles_V{iSubj, iToPlot, Sorting_Cdt, iROI};
                X_sort(iSubj, :) = All_X_sort_V{iSubj, iToPlot, Sorting_Cdt, iROI};
                Profiles(:, :, iSubj) = All_Profiles_V{iSubj, iToPlot, iCdt, iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            % computes vertex wise correlation / regression between predictor
            % and predicted rasters plots
            [rho,slope] = CorRegRaster(Profiles, DesMat, iToPlot, X_sort);
            
            
            %% plot main raster
            subplot(1, NbSubPlot, 2:3)
            PlotRectangle(NbLayers, 10, 1, 0)
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            subplot(1, NbSubPlot, 2:3)
            hold on
            colormap(ColorMap);
            
            % plot raster (add some gaussian smoothing across vertices)
            raster_to_plot = mean(imgaussfilt(Profiles,[size(Profiles,1)/100 .0001]), 3);
            if demean_raster && strcmp(ToPlot{iToPlot}, 'Lin')
                mnval = mean(raster_to_plot, 2);
                for i=1:size(raster_to_plot, 1)
                    raster_to_plot(i,:) = raster_to_plot(i,:) - mnval(i);
                end
            end
            imagesc(raster_to_plot, CLIM)
            axis([0.5 NbLayers+.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            
            %% plot sorting raster
            subplot(1, NbSubPlot, 1)
            PlotRectangle(NbLayers, 10, 1, 0)
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            subplot(1, NbSubPlot, 1)
            hold on
            colormap(ColorMap);
            
            % plot raster (add some gaussian smoothing across vertices)
            raster_to_plot = mean(imgaussfilt(Sorting_Raster,[size(Sorting_Raster,1)/100 .0001]),3);
            if demean_raster && strcmp(ToPlot{iToPlot}, 'Lin')
                mnval = mean(raster_to_plot, 2);
                for i=1:size(raster_to_plot, 1)
                    raster_to_plot(i,:) = raster_to_plot(i,:) - mnval(i);
                end
            end
            imagesc(raster_to_plot, CLIM2)
            axis([0.5 NbLayers+.5 0 size(Profiles,1)])
            
            set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                'ytick', [],'yticklabel', [])
            
            
            %% plot sorting value
            ax = gca;
            YLabel = [];
            PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, 0)
            set(gca,'fontsize', 20)
            
            % plot regression coefficient
            if plot_reg_coeff
                subplot(1, NbSubPlot, 4) %#ok<*UNRCH>
                ax = gca;
                axis off
                PlotCorrCoeff(ax, slope, ToPlot{iToPlot}, 0, 0, ax.Position(3), ax.Position(4), ...
                    [0.9 1.6 -0.1 0.5], ToPermute)
                set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], 'fontsize',10,'yaxislocation', 'right')
            end
            
            
            %% saves table results and figure
            SavedTxt = fullfile(FigureFolder, 'Baseline', [fig_name '.csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, ...
                [Name{iCdt} '_fV'], slope, ToPermute)
            
            
            print(gcf, fullfile(FigureFolder,'Baseline', [fig_name '.tif']), '-dtiff')
            
        end
        
    end
    
end
