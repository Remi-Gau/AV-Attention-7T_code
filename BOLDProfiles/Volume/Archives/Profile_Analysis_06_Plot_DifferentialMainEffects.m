clear; clc; close all

NbLayers = 6;

StartFolder=pwd;
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '14';...
    '15';...
    '16'
    ];

FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles', strcat(num2str(NbLayers), '_layers'));
cd(FigureFolder)
load(strcat('Data_Block_', num2str(NbLayers), '_Layers', '.mat'), 'AllSubjects_Data')


%% Plots
for ROI_Ind = [7 10 22] %1:length(AllSubjects_Data)
    
    close all
    
    cd(FigureFolder)
    mkdir(AllSubjects_Data(ROI_Ind).name)
    cd(fullfile(FigureFolder, AllSubjects_Data(ROI_Ind).name))
    
    fprintf([AllSubjects_Data(ROI_Ind).name '\n'])
    
    Include = find(AllSubjects_Data(ROI_Ind).Include);
%     Include = ~any(AllSubjects_Data(ROI_Ind).VoxelCount(:,2:end)<300,2);
    
    if ~isempty(Include)
        

        %% Plot Differentials Main effect +/- SEM with subjects
        Visible='off';
        PRINT=4;
        
        Name = [AllSubjects_Data(ROI_Ind).name '-Block-Differentials-main-effects-with-subjects'];
        ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Differential.MainEffect.MEAN);
        ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Differential.MainEffect.SEM);
        ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.MainEffect.DATA;
        ToPlot.MAX = max(ToPlot.Subjects.Data(:));
        ToPlot.MIN = min(ToPlot.Subjects.Data(:));
        ToPlot.Subjects.List = SubjectList(Include,:);
        ToPlot.Subjects.Legend.Bool = 0;
        ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
        ToPlot.Include=Include;
        
        ToPlot.Subjects.SubPlots = 1:2;
        
        ToPlot.P = AllSubjects_Data(ROI_Ind).Differential.Blocks.MainEffects.Beta.P;
        
        ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
        
        PlotLayers(ToPlot, Name, Visible, PRINT)
        
        clear ToPlot
        
        
        %% Plot Differentials restricts Main effect +/- SEM with subjects
        Visible='off';
        PRINT=4;
        
        Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Differentials-main-effects-with-subjects'];
        ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.MainEffect.MEAN);
        ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.MainEffect.SEM);
        ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.MainEffect.DATA;
        ToPlot.MAX = max(ToPlot.Subjects.Data(:));
        ToPlot.MIN = min(ToPlot.Subjects.Data(:));
        ToPlot.Subjects.List = SubjectList(Include,:);
        ToPlot.Subjects.Legend.Bool = 0;
        ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
        ToPlot.Include=Include;
        
        ToPlot.Subjects.SubPlots = 1:2;
        
        ToPlot.P = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.MainEffects.Beta.P;
        
        ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
        
        PlotLayers(ToPlot, Name, Visible, PRINT)
        
        clear ToPlot

    end
    
end

cd(StartFolder)