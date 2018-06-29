clear; clc;

NbLayers = 6;

StartFolder=fullfile(pwd, '..','..');
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


ToPlot.shadedErrorBar.Legend.Legend = {};
ToPlot.shadedErrorBar.mean = [];
ToPlot.shadedErrorBar.errorbar = [];

j=0;
 
%% Plots
for ROI_Ind = [11:13 20:21 7 10] %[22 20 21 7 10] %[11:13 20:21] %1:length(AllSubjects_Data)
    
    close all
    
    j = j+1;
    
    fprintf([AllSubjects_Data(ROI_Ind).name '\n'])
    
    Include = find(AllSubjects_Data(ROI_Ind).Include);
    
    Name = 'ProfileAcrossRegions2';
    
    ToPlot.shadedErrorBar.mean(:,:,j) = flipud(AllSubjects_Data(ROI_Ind).Differential.MEAN);
    ToPlot.shadedErrorBar.errorbar(:,:,j) = flipud(AllSubjects_Data(ROI_Ind).Differential.SEM);
    
    ToPlot.shadedErrorBar.Legend.Legend{j} = AllSubjects_Data(ROI_Ind).name;
    
    
end

%%
Visible='on';
PRINT=1;

ToPlot.shadedErrorBar.Legend.Bool=1;
ToPlot.Include=Include;

tmp = ToPlot.shadedErrorBar.mean+ToPlot.shadedErrorBar.errorbar;
ToPlot.MAX = max(tmp(:));
ToPlot.Max = max(max(max(tmp(:,[4 8:12],:))));
tmp = ToPlot.shadedErrorBar.mean-ToPlot.shadedErrorBar.errorbar;
ToPlot.MIN = min(tmp(:));
ToPlot.Min = min(min(min(tmp(:,[4 8:12],:))));

PlotLayers(ToPlot, Name, Visible, PRINT)


cd(StartFolder)
