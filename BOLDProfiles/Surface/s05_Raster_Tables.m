% plot rasters (and corresponding tables) used for the paper
% will actually plot some extra

close all
clear
clc

StartFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T';
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox')

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

load(fullfile(StartFolder,'Results','Profiles','Surfaces','Raster_CV.mat'))

Subj2Include = true(11,1);

for iSubj=1:sum(Subj2Include)
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];

ToPlot={'Cst','Lin'};

NbLayers = 6;
DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


%% Grp level  ; raster Att = f(Percentile of V)
mkdir(fullfile(FigureFolder,'Attention'));
NameSwitch={'Att_V-A', 'Astim_Att_V-A', 'Vstim_Att_V-A'};

for iToPlot = 1:2

    for iCdt = 1:3
        
        for iROI = 1:4

            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_A,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_AttfV{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_AttfV{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);

            SavedTxt = fullfile(FigureFolder,'Attention', ...
                ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '.csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, NameSwitch{iCdt}, slope, ToPermute)
            
        end
        
    end

end


%% Grp level  ; raster CrossMod = f(Percentile of V stim)
mkdir(fullfile(FigureFolder,'CrossSensory'));
Name = {'AV-A';'AV-V'};

for iToPlot = 1:2
    
    for iCdt = 1:2
        
        for iROI = 1:4

            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_V{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_Cross_V{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Cross_V{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);

            SavedTxt = fullfile(FigureFolder,'CrossSensory', ...
                ['GrpLvl_Raster_' Name{iCdt} '_fV_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, [Name{iCdt} '_fV'], slope, ToPermute)

        end
        
    end

end


%% Grp level  ; raster Att = f(Percentile of A)
mkdir(fullfile(FigureFolder,'Attention'));
% Name={'Att_A-V', 'Astim_Att_A-V', 'Vstim_Att_A-V'};
NameSwitch={'Att_V-A', 'Astim_Att_V-A', 'Vstim_Att_V-A'};

for iToPlot = 1:2
    
    for iCdt = 1:3
        
        for iROI = 1:4

            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_Profiles_A,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,1,iROI};
                X_sort(iSubj,:) = All_X_sort_AttfA{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_AttfA{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);
            
            SavedTxt = fullfile(FigureFolder,'Attention', ...
                ['Raster_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '_' NameSwitch{iCdt} '_f(A).csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, NameSwitch{iCdt}, slope, ToPermute)

        end
        
    end

end


%% Grp level  ; raster CrossMod = f(Percentile of V stim)
close all
mkdir(fullfile(FigureFolder,'CrossSensory'));
Name = {'AV-A';'AV-V'};

for iToPlot = 1:2
    
    for iCdt = 1:2
        
        for iROI = 1:4

            clear X_sort Profiles Sorting_Raster
            
            for iSubj = 1:size(All_X_sort_Cross,1)
                Sorting_Raster(:,:,iSubj) = All_Profiles_A{iSubj,iToPlot,2,iROI};
                X_sort(iSubj,:) = All_X_sort_Cross_A{iSubj,iToPlot,iCdt,iROI};
                Profiles(:,:,iSubj) = All_Profiles_Cross_A{iSubj,iToPlot,iCdt,iROI};
            end
            
            X_sort = X_sort(Subj2Include,:);
            Profiles = Profiles(:,:,Subj2Include);
            
            [rho,slope]=CorRegRaster(Profiles,DesMat,iToPlot,X_sort);

            SavedTxt = fullfile(FigureFolder,'CrossSensory', ...
                ['GrpLvl_Raster_' Name{iCdt} '_fA_Final_' ToPlot{iToPlot} '_CV_' ROI(iROI).name '.csv']);
            PrintTableCorrCoeff(SavedTxt, ROI(iROI).name, ToPlot{iToPlot}, [Name{iCdt} '_fA'], slope, ToPermute)

        end
        
    end
    
end