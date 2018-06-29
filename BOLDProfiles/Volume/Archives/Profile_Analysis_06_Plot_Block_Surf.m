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
load(strcat('Data_Block_Surf_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat'), 'AllSubjects_Data')


%% Plots
for ROI_Ind = 1:length(AllSubjects_Data)
    
    close all
    
    cd(FigureFolder)
    mkdir(AllSubjects_Data(ROI_Ind).name)
    cd(fullfile(FigureFolder, AllSubjects_Data(ROI_Ind).name))
    
    fprintf([AllSubjects_Data(ROI_Ind).name '\n'])
    
    Include = find(AllSubjects_Data(ROI_Ind).Include);
%     Include = find(all([~any(AllSubjects_Data(ROI_Ind).VoxelCount(:,2:end)<500,2) ...
%         AllSubjects_Data(ROI_Ind).Include],2));
    
    if ~isempty(Include)
        
        %% Plot Mean +/- SEM
                Visible='on';
                PRINT=3;
        
                Name = [AllSubjects_Data(ROI_Ind).name '-Surf-Block'];
                ToPlot.MAX = max(AllSubjects_Data(ROI_Ind).MEAN(:) + ...
                    AllSubjects_Data(ROI_Ind).SEM(:));
                ToPlot.MIN = min(AllSubjects_Data(ROI_Ind).MEAN(:)-...
                    AllSubjects_Data(ROI_Ind).SEM(:));
                
                
                ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MEAN);
                ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).SEM);
        
                PlotLayers(ToPlot, Name, Visible, PRINT)
        
                clear ToPlot
        
        %% Plot Mean +/- SEM with subjects
                Visible='on';
                PRINT=3;
        
                Name = [AllSubjects_Data(ROI_Ind).name '-Surf-Block-With subjects'];
                ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MEAN);
                ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).SEM);
                ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DATA(:,:,Include);
                ToPlot.MAX = max(ToPlot.Subjects.Data(:));
                ToPlot.MIN = min(ToPlot.Subjects.Data(:));
                ToPlot.Subjects.List = SubjectList(Include,:);
                ToPlot.Subjects.Legend.Bool = 0;
                ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
                ToPlot.Subjects.SubPlots = 1:6;
                ToPlot.Include=Include;
        
                PlotLayers(ToPlot, Name, Visible, PRINT)
        
                clear ToPlot
        
        %% Plot differential +/- SEM
%                 Visible='off';
%                 PRINT=1;
%         
%                 Name = [AllSubjects_Data(ROI_Ind).name '-Block-Differential'];
%                 ToPlot.MAX = max(AllSubjects_Data(ROI_Ind).Differential.MEAN(:) + ...
%                     AllSubjects_Data(ROI_Ind).Differential.SEM(:));
%                 ToPlot.MIN = min(AllSubjects_Data(ROI_Ind).Differential.MEAN(:)-...
%                     AllSubjects_Data(ROI_Ind).Differential.SEM(:));
%                 ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Differential.MEAN);
%                 ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Differential.SEM);
%         
%                 tmp = (AllSubjects_Data(ROI_Ind).Differential.MEAN + ...
%                     AllSubjects_Data(ROI_Ind).Differential.SEM);
%                 tmp = tmp(:,[4 8 9:12]);
%                 ToPlot.Max = max(tmp(:));
%         
%                 tmp = (AllSubjects_Data(ROI_Ind).Differential.MEAN-...
%                     AllSubjects_Data(ROI_Ind).Differential.SEM);
%                 tmp = tmp(:,[4 8 9:12]);
%                 ToPlot.Min = min(tmp(:));
%         
%                 ToPlot.Include=Include;
%         
%                 ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%                 PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%                 clear ToPlot
        
        
        %% Plot Derivative +/- SEM
%                 Visible='off';
%                 PRINT=1;
%         
%                 Name = [AllSubjects_Data(ROI_Ind).name '-Block-Derivative'];
%                 ToPlot.MAX = max(AllSubjects_Data(ROI_Ind).Derivative.MEAN(:) + ...
%                     AllSubjects_Data(ROI_Ind).Derivative.SEM(:));
%                 ToPlot.MIN = min(AllSubjects_Data(ROI_Ind).Derivative.MEAN(:)-...
%                     AllSubjects_Data(ROI_Ind).Derivative.SEM(:));
%                 ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Derivative.MEAN);
%                 ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Derivative.SEM);
%         
%                 PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%                 clear ToPlot
        
         %% Plot Derivative +/- SEM with subjects
%                 Visible='off';
%                 PRINT=1;
%         
%                 Name = [AllSubjects_Data(ROI_Ind).name '-Block-Derivative with subjects'];
%                 ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Derivative.MEAN);
%                 ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Derivative.SEM);
%                 ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Derivative.DATA;
%                 ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%                 ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%                 ToPlot.Subjects.List = SubjectList(Include,:);
%                 ToPlot.Subjects.Legend.Bool = 0;
%                 ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%                 ToPlot.Subjects.SubPlots = 1:12;
%                 ToPlot.Include=Include;
%         
%                 PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%                 clear ToPlot
        
        %% Plot differential +/- SEM with subjects
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Differential with subjects'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Differential.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Differential.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).Differential.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot differential +/- SEM with subjects (but not for differential plots)
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Differential with subjects 2'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Differential.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Differential.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = [1:3 5:7];
%         ToPlot.Include=Include;
%         
%         tmp = ToPlot.shadedErrorBar.mean+ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Max = max(tmp(:));
%         tmp = ToPlot.shadedErrorBar.mean-ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot differential +/- SEM with subjects with legend
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Differential with subjects-Legend'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).Differential.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).Differential.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).Differential.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        %% Betas
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Betas'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Betas with legend
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Betas-Legend'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        
        
        %% Plot restricted differential +/- SEM with subjects
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Differential with subjects'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).DifferentialRestrict.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot restricted differential +/- SEM with subjects (but not for differential plots)
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Differential with subjects 2'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = [1:3 5:7];
%         ToPlot.Include=Include;
%         
%         tmp = ToPlot.shadedErrorBar.mean+ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Max = max(tmp(:));
%         tmp = ToPlot.shadedErrorBar.mean-ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot restricted differential +/- SEM with subjects with legend
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Differential with subjects-Legend'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).DifferentialRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).DifferentialRestrict.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        %% Betas restricted
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Betas'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Betas restricted with legend
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Betas-Legend'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DifferentialRestrict.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        

        %% Plot Main effect +/- SEM with subjects
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Main effect with subjects'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffects.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffects.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffects.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).MainEffects.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot Main effect +/- SEM with subjects (but not for differential plots)
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Main effect with subjects 2'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffects.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffects.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffects.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = [1:3 5:7];
%         ToPlot.Include=Include;
%         
%         tmp = ToPlot.shadedErrorBar.mean+ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Max = max(tmp(:));
%         tmp = ToPlot.shadedErrorBar.mean-ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot Main effect +/- SEM with subjects with legend
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Main effect with subjects-Legend'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffects.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffects.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffects.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).MainEffects.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        %% Betas Main effect
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Betas Main effect'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Betas Main effect with legend
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-Block-Betas Main effect-Legend'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffects.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        
        
        %% Plot restricted Main effect +/- SEM with subjects
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Main effect with subjects'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).MainEffectsRestrict.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot restricted Main effect +/- SEM with subjects (but not for differential plots)
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Main effect with subjects 2'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = [1:3 5:7];
%         ToPlot.Include=Include;
%         
%         tmp = ToPlot.shadedErrorBar.mean+ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Max = max(tmp(:));
%         tmp = ToPlot.shadedErrorBar.mean-ToPlot.shadedErrorBar.errorbar;
%         tmp = tmp(:,[4 8 9:11],:);
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.P;
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Plot restricted Main effect +/- SEM with subjects with legend
%         Visible='off';
%         PRINT=4;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Main effect with subjects-Legend'];
%         ToPlot.shadedErrorBar.mean = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.MEAN);
%         ToPlot.shadedErrorBar.errorbar = flipud(AllSubjects_Data(ROI_Ind).MainEffectsRestrict.SEM);
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.DATA;
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         tmp = (AllSubjects_Data(ROI_Ind).MainEffectsRestrict.DATA);
%         tmp = tmp(:,[4 8 9:12],:);
%         ToPlot.Max = max(tmp(:));
%         ToPlot.Min = min(tmp(:));
%         
%         ToPlot.VoxelCount = AllSubjects_Data(ROI_Ind).VoxelCount(Include,2:end);
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotLayers(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        
        %% Betas Main effect restricted
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Betas Main effect'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 0;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot
        
        %% Betas Main effect restricted with legend
%         Visible='off';
%         PRINT=1;
%         
%         Name = [AllSubjects_Data(ROI_Ind).name '-BlockRestricted-Betas Main effect-Legend'];
%         
%         ToPlot.MEAN = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.MEAN;
%         ToPlot.SEM = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.SEM;
%         ToPlot.P = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.P;
%         ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).MainEffectsRestrict.Blocks.Beta.DATA(:,:,Include);
%         
%         
%         ToPlot.MAX = max(ToPlot.Subjects.Data(:));
%         ToPlot.MIN = min(ToPlot.Subjects.Data(:));
%         
%         ToPlot.Subjects.List = SubjectList(Include,:);
%         ToPlot.Subjects.Legend.Bool = 1;
%         ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
%         ToPlot.Subjects.SubPlots = 1:12;
%         ToPlot.Include=Include;
%         
%         ToPlot.MainEffect = 1;
%         
%         PlotBetas(ToPlot, Name, Visible, PRINT)
%         
%         clear ToPlot



    end
    
end

cd(StartFolder)


%%
% FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles');
%
% cd(FigureFolder)
%
% List = dir('*.pdf');
%
% Command = [];
% for iFile = 1:numel(List)
%     Command = [Command ' ' List(iFile).name]; %#ok<AGROW>
% end
%
%
% system([...
%     'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
%     '-sOutputFile=' fullfile(FigureFolder, 'ProfilesVolumeAllSubjects.pdf') ' ' Command])
