%%
clc
close all
clear

StartFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives';
addpath(genpath(fullfile(StartFolder, 'SubFun')))

Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

% Color map
X = 0:0.001:1;
ColorMap = seismic(size(X,2));
clear R G B

FigDim = [100, 100, 1000, 1500];
Visibility = 'on';

% StartFolder='/media/rxg243/BackUp2/AV_Integration_7T_2/Scripts/BOLDProfiles/Surface/../../..';
% addpath(genpath(fullfile(StartFolder, 'SubFun')))

StartFolder=fullfile(pwd, '..','..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

mkdir(fullfile(FigureFolder,'Baseline', 'RasterAllCdt'))

% load(fullfile(StartFolder,'Results','Profiles','Surfaces','RasterAllCdt_Lin_CV.mat'))
% load(fullfile(StartFolder,'Results','Profiles','Surfaces','RasterAllCdt_CV.mat'))
load(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster_CV.mat'))

Color = [...
    0,0,1;
    .5,.5,1;
    1,.5,.5;
    1,0,0];

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};





load(fullfile(StartFolder,'Figures','8_layers','MinNbVert.mat'),'MinVert')


%% Grp level  ; raster stim = f(Percentile of V)
close all
CLIM = [-1 1];
nMax = 30;
% Cdt =[2 1;2 2;2 3;2 4;2 5;2 6];
Cdt =[2 2;2 1];
Name={'A','V','AV','A','V','AV'};
ToPlot={'Cst','Lin'};

for iToPlot = 1%:numel(ToPlot)
    
    close all
    
    for iROI = 2 % numel(ROI)
        
        NbBin = MinVert(strcmp(ROI(iROI).name,{MinVert.name}')).MinVert;
        
        Splits = floor(linspace(1,NbBin,nMax+1));
        
        fprintf('    %s\n',ROI(iROI).name)
        
        for n=1:nMax
            
            FileName = fullfile(FigureFolder,'Baseline', 'RasterAllCdt', ...
                ['GrpLvl_Raster_CrossMod_f_V_stim_' ToPlot{iToPlot} '_' ROI(iROI).name '_CV.gif']);
            
            h = figure('name', ToPlot{iToPlot}, 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility); %#ok<*UNRCH>
            
            for iCdt = 1:2
                
                clear X_sort Profiles Sorting_Raster
                
                for iSubj = 1:size(All_Profiles_V,1)
                    Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                    X_sort(iSubj,:) = All_X_sort_V{iSubj,iToPlot,Cdt(iCdt,2),iROI};
                    if iCdt==2
                        Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,Cdt(iCdt,2),iROI};
                    else
                        Profiles(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,Cdt(iCdt,2),iROI};
                    end
                end
                
%                 Sorting_Raster = Sorting_Raster(:,:,Subj2Include);
                
%                 X_sort = X_sort(Subj2Include,:);
%                 Profiles = Profiles(:,:,Subj2Include);
                MeanProfiles = mean(Profiles,3);
                
                
                subplot(1,2,iCdt)
                
                colormap(ColorMap);
                imagesc(flipud(MeanProfiles), CLIM)
                imagesc(imgaussfilt(flipud(MeanProfiles),[size(MeanProfiles,1)/200 .001]), CLIM)
                %             imagesc(MeanProfiles, [-5 5])
                
                axis([0.5 6.5 0 size(Profiles,1)])
                
                set(gca,'tickdir', 'out', 'xtick', [],'xticklabel', [], ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                if iCdt==1
                    t = title('[V - Fix]');
                else
                    t = title('[AV - V]');
                end
                
                set(t,'fontsize',10)
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'ytick', [],'yticklabel', [], ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                
                
                ax = gca;
                axes('Position',ax.Position);
                
                hold on
                
                errorbar(1:6,mean(MeanProfiles(Splits(n):Splits(n+1),:)),...
                    nansem(squeeze(mean(Profiles(Splits(n):Splits(n+1),:,:))),2), 'k', 'linewidth', 2)
                for iSubj=1:size(Profiles,3)
                    plot(1:6,mean(Profiles(Splits(n):Splits(n+1),:,iSubj)), ':k', 'linewidth', .5)
                end
                plot([1 6],[0 0], '--k')
                
                axis([0.5 6.5 -2 2])
                
                DephLevels = round(linspace(100,0,8));
                DephLevels([1;end]) = [];
                set(gca,'color', 'none', 'tickdir', 'out', 'xtick', 1:6,'xticklabel',  DephLevels, ...
                    'YAxisLocation','right', 'ytick', -10:10,'yticklabel', -10:10, ...
                    'ticklength', [0.01 0], 'fontsize', 10)
                
                t=xlabel('cortical depth');
                set(t,'fontsize',10)
                
                
                if iCdt==1
                    ax = gca;
                    
                    YLabel = sprintf('%s\nPerc %s %s stim',...
                        strrep(ROI(iROI).name,'_','-'),ToPlot{iToPlot}, Name{Cdt(iCdt,1)});
%                     PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, Sorting_Raster, CLIM)
                      PlotSortedValues(ax, X_sort, NbBin, Profiles, YLabel, 1, [], CLIM, 0)
                  
                    
                    tmp = axis;
                    r=rectangle('Position', [tmp(1) Splits(n) tmp(2)-tmp(1) (Splits(n+1)-Splits(n))]);
                    set(r,'linewidth',2, 'EdgeColor', 'r')
                    
                    clear tmp
                    
                end
                
                if iCdt==3
                    PlotColorBar(ax, ColorMap, CLIM)
                end
                
            end
            
            mtit([strrep(ROI(iROI).name,'_',' ') ' - Percentile V stim - ' ToPlot{iToPlot}], 'fontsize', 14, 'xoff',0,'yoff',.025)
            
            pause(0.2)
            
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            if n == 1
                imwrite(imind,cm,FileName,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,FileName,'gif','WriteMode','append','DelayTime',1);
            end
            
            close all
            
        end
        
        close all
        
    end
    
end

