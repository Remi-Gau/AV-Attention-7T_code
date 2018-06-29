clear; clc; close all

NbLayers = 3;

StartFolder=fullfile(pwd,'..','..');
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
for ROI_Ind = 1:length(AllSubjects_Data)
    
    close all
    
    cd(FigureFolder)
    mkdir(AllSubjects_Data(ROI_Ind).name)
%     cd(fullfile(FigureFolder, AllSubjects_Data(ROI_Ind).name))
    
    fprintf([AllSubjects_Data(ROI_Ind).name '\n'])
    
    Include = find(AllSubjects_Data(ROI_Ind).Include);
%     Include = ~any(AllSubjects_Data(ROI_Ind).VoxelCount(:,2:end)<300,2);
    
    if ~isempty(Include)
        
        %% Plot Bi vs Uni +/- SEM with subjects with legend
        Visible='off';
        PRINT=4;
        
        Name = [AllSubjects_Data(ROI_Ind).name '-Block-BiVSUniOnly-' num2str(NbLayers) '_Layers'];
        ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).BiVSUni.MEAN(:,[7:8]));
        ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).BiVSUni.SEM(:,[7:8]));
        ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).BiVSUni.DATA(:,[7:8],:);
        ToPlot.MAX = max(ToPlot.Subjects.Data(:));
        ToPlot.MIN = min(ToPlot.Subjects.Data(:));
        ToPlot.Subjects.List = SubjectList(Include,:);
        ToPlot.Subjects.Legend.Bool = 0;
        ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
        ToPlot.Subjects.SubPlots = 1:2;
        ToPlot.Include=Include;
        
        ToPlot.P = AllSubjects_Data(ROI_Ind).BiVSUni.Blocks.Beta.P(:,[7:8]);
        
      
        PlotLayers(ToPlot, Name, Visible, PRINT)
        
        clear ToPlot
        
    end
    
end

cd(StartFolder)