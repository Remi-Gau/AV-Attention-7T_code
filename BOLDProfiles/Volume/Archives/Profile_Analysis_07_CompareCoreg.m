clear; close all; clc;

StartFolder=pwd;

addpath(genpath(fullfile(StartFolder,'SubFun')))

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
    %     '15';...
    '16'
    ];

MinLayer = 1;
NbLayers = 10;
TargetLayerFile = ['T1_' sprintf('%02.0f', NbLayers) '_Layers'];

Include = logical([1 1 1 1 0 1 1 1 1 1 1 1]);

AllSubjects_Data = struct();

for SubjInd=1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    ROI_Folder = fullfile(SubjectFolder, 'ROI_MNI');
    
    AnalysisFolder = fullfile(SubjectFolder, 'Transfer', 'ROI');
    [~,~,~]=mkdir(AnalysisFolder);
    
    %% Copy ROIs
    cd(AnalysisFolder)
    
    copyfile(fullfile(ROI_Folder,'rwROI_auditory_*_MNI.nii'), fullfile(AnalysisFolder));
    
    movefile('rwROI_auditory_Te10_MNI.nii', 'TE_1.0_MNI.nii');
    movefile('rwROI_auditory_Te11_MNI.nii', 'TE_1.1_MNI.nii');
    movefile('rwROI_auditory_Te12_MNI.nii', 'TE_1.2_MNI.nii');
    
    copyfile(fullfile(ROI_Folder,'rwROI_Visual_*_MNI.nii'), fullfile(AnalysisFolder));
    
    movefile('rwROI_Visual_hOC1_MNI.nii', 'V1_MNI.nii')
    movefile('rwROI_Visual_hOc2_MNI.nii', 'V2_MNI.nii')
    movefile('rwROI_Visual_hOc3d_MNI.nii', 'V3d_MNI.nii')
    movefile('rwROI_Visual_hOc3v_MNI.nii', 'V3v_MNI.nii')
    movefile('rwROI_Visual_hOc4d_MNI.nii', 'V4d_MNI.nii')
    movefile('rwROI_Visual_hOc4v_MNI.nii', 'V4v_MNI.nii')
    movefile('rwROI_Visual_hOc5_MNI.nii', 'V5_MNI.nii')
    
    delete('rw*')
    
    
    Files2Reorient = {};
    tmp=dir('*.nii');
    for iROI=1:numel(tmp)
        Files2Reorient{end+1,1} = fullfile(AnalysisFolder,tmp(iROI).name);
    end
    clear tmp i
    
    
    %% Reorient ROIs
    cd(fullfile(SubjectFolder, 'Transfer'));
    
    MatFiles = dir('*mat');
    MatFiles = char(MatFiles.name);
    
    Reorient = strcmp(cellstr(MatFiles(:,1:14)),'ReorientMatrix');
    if any(Reorient)
        if sum(Reorient)>1
            error('Images have been reoriented more than once.')
            return
        else
            load(deblank(MatFiles(Reorient,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    clear Reorient M
    
    Coregister = strcmp(cellstr(MatFiles(:,1:11)),'CoregMatrix');
    if any(Coregister)
        if sum(Coregister)>1
            error('Images have been coregistered more than once.')
            return
        else
            load(deblank(MatFiles(Coregister,:)), 'M')
            spm_reorient(Files2Reorient, M)
        end
    end
    
    clear Files2Reorient M
    
    
    %% Reslice ROIs
    cd(AnalysisFolder)
    if isempty(dir('r*.nii'))
        
        flags = struct(...
            'mask', 0, ...
            'mean', 0, ...
            'interp', 1, ...
            'which', 1, ...
            'wrap', [0 0 0], ...
            'prefix', 'r' ...
            );
        
        ImagesFiles2Process = {};
        
        ImagesFiles2Process{1,1} = fullfile(SubjectFolder, 'Transfer', 'rcon_0001.nii');
        
        tmp=dir('*.nii');
        for iROI=1:numel(tmp)
            ImagesFiles2Process{end+1,1} = fullfile(AnalysisFolder,tmp(iROI).name); %#ok<SAGROW>
        end
        clear tmp i
        
        spm_reslice(ImagesFiles2Process,flags)
        
    end
    
    
    %% Average
    cd(AnalysisFolder)
    ROIs=dir('r*.nii');
    
    if Include(SubjInd)
        
        if exist(fullfile(AnalysisFolder, ['Data_' num2str(NbLayers) '_Layers' '.mat']), 'file')
            load(fullfile(AnalysisFolder, ['Data_' num2str(NbLayers) '_Layers' '.mat']), 'DATA')
            
            for iROI=1:numel(ROIs)
                disp(ROIs(iROI).name)
                AllSubjects_Data(iROI).name=DATA(iROI).name;
                AllSubjects_Data(iROI).DATA(:,:,SubjInd) = DATA(iROI).DATA;
                AllSubjects_Data(iROI).DATA2(:,:,SubjInd) = DATA(iROI).DATA2;
            end
        else
            
            LayersVolumeStruc = spm_vol(fullfile(SubjectFolder, 'Structural', 'CBS', ...
                'Layering',[TargetLayerFile '.nii']));
            
            cd(fullfile(SubjectFolder, 'Transfer'))
            
            ListOfContrastFiles = dir('rcon*def.nii');
            SourceVolumeStruc = char(ListOfContrastFiles.name);
            SourceVolumeStruc = spm_vol(SourceVolumeStruc);
            
            ListOfContrastFiles = dir('rcon*.nii');
            SourceVolumeStruc2={};
            for i=1:numel(ListOfContrastFiles)
                if ListOfContrastFiles(i).name(end-4)~='f'
                    SourceVolumeStruc2{end+1,1} = ListOfContrastFiles(i).name; %#ok<SAGROW>
                end
            end
            SourceVolumeStruc2 = char(SourceVolumeStruc2);
            SourceVolumeStruc2 = spm_vol(SourceVolumeStruc2);
            
            for iROI=1:numel(ROIs)
                
                disp(ROIs(iROI).name)
                
                AllSubjects_Data(iROI).name = ROIs(iROI).name;
                
                DATA(iROI).name = ROIs(iROI).name; %#ok<*SAGROW>
                
                ROI = spm_read_vols(spm_vol(fullfile(AnalysisFolder,...
                    ROIs(iROI).name)));
                
                [X, Y, Z] = ind2sub(size(ROI), find(ROI));
                XYZ = [X'; Y'; Z'];
                clear X Y Z
                
                Vox = spm_get_data(LayersVolumeStruc,XYZ);
                
                cd(fullfile(SubjectFolder, 'Transfer'))
                Vol = spm_get_data(SourceVolumeStruc,XYZ);
                Vol2 = spm_get_data(SourceVolumeStruc2,XYZ);
                
                
                for iLayer=1:NbLayers
                    for iCdti=1:6
                        AllSubjects_Data(iROI).DATA(iLayer,iCdti,SubjInd) = nanmean(Vol(iCdti,Vox==iLayer));
                        AllSubjects_Data(iROI).DATA2(iLayer,iCdti,SubjInd) = nanmean(Vol2(iCdti,Vox==iLayer));
                        DATA(iROI).DATA(iLayer,iCdti) = nanmean(Vol(iCdti,Vox==iLayer));
                        DATA(iROI).DATA2(iLayer,iCdti) = nanmean(Vol2(iCdti,Vox==iLayer));
                    end
                end
                
            end
            save(fullfile(AnalysisFolder, ['Data_' num2str(NbLayers) '_Layers' '.mat']), 'DATA')
        end
    end
    
end

%% Plot Mean +/- SEM with subjects
cd(StartFolder)
delete('*.tif')

ROIs={...
    'TE_1.0_MNI';...
    'TE_1.1_MNI';...
    'TE_1.2_MNI';...
    'V1_MNI';...
    'V2_MNI';...
    'V3d_MNI';...
    'V3v_MNI';...
    'V4d_MNI';...
    'V4v_MNI';...
    'V5_MNI'};

Visible = 'off';
PRINT = '100';

Command = [];

for iROI=1:numel(ROIs)
    
    %%
    disp(ROIs{iROI})
    
    Idx = find(strcmp({AllSubjects_Data.name},['r' ROIs{iROI} '.nii']));
        
    DATA = AllSubjects_Data(Idx).DATA(:,:,Include);
    
    MEAN = nanmean(DATA,3);
    SEM = nansem(DATA,3);
    STD = nanstd(DATA,3);
    
    A = DATA;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential.DATA = A;
    Differential.MEAN = nanmean(A,3);
    Differential.STD = nanstd(A,3);
    Differential.SEM = nansem(A,3);
    clear A
    
    
    
    DATA2 = AllSubjects_Data(Idx).DATA2(:,:,Include);
    
    MEAN2 = nanmean(DATA2,3);
    SEM2 = nansem(DATA2,3);
    STD2 = nanstd(DATA2,3);
    
    A = DATA2;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential2.DATA = A;
    Differential2.MEAN = nanmean(A,3);
    Differential2.STD = nanstd(A,3);
    Differential2.SEM = nansem(A,3);
    clear A
    
    
    MAX = max([DATA(:);DATA2(:)]);
    MIN = min([DATA(:);DATA2(:)]);
    
    tmp = Differential.DATA;
    tmp = tmp(:,[4 8 9:12],:);
    tmp2 = Differential2.DATA;
    tmp2 = tmp2(:,[4 8 9:12],:);
    Max = max([tmp(:);tmp2(:)]);
    Min = min([tmp(:);tmp2(:)]);
    
    %%
    Name = [ROIs{iROI} ' SyN - With subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM );
    ToPlot.Subjects.Data = DATA;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
%     PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
    
    %%
    Name = [ROIs{iROI} ' SPM - With subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN2);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM2);
    ToPlot.Subjects.Data = DATA2;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
%     PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
    
    %%    
    Name = [ROIs{iROI} ' SyN - Differential with subjects'];
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.shadedErrorBar.mean = flipud(Differential.MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(Differential.SEM);
    ToPlot.Subjects.Data = Differential.DATA;
    
    
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    ToPlot.Max = Max;
    ToPlot.Min = Min;
    
    ToPlot.Include=Include;

    PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
    
    %%
    Name = [ROIs{iROI} ' SPM - Differential with subjects'];
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.shadedErrorBar.mean = flipud(Differential2.MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(Differential2.SEM);
    ToPlot.Subjects.Data = Differential2.DATA;
    
    
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
    ToPlot.Max = Max;
    ToPlot.Min = Min;
    
    ToPlot.Include=Include;

    PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
    
end

AllSubjects_Data_Ori = AllSubjects_Data;

%%
FigureFolder = fullfile(StartFolder, 'ROI_Analysis', 'Figures', strcat(num2str(NbLayers), '_layers'));
cd(FigureFolder)
load(strcat('Data_', num2str(NbLayers), '_Layers', '.mat'), 'AllSubjects_Data')

cd(StartFolder)

for iROI = 3:5
    disp(AllSubjects_Data(5).SubROI(iROI).name)
    
    Idx = find(strcmp({AllSubjects_Data_Ori.name},...
        ['r' AllSubjects_Data(5).SubROI(iROI).name '_MNI.nii']));
    DATA = AllSubjects_Data_Ori(Idx).DATA(:,:,Include);
    DATA2 = AllSubjects_Data_Ori(Idx).DATA2(:,:,Include);
    MAX = max([DATA(:);DATA2(:)]);
    MIN = min([DATA(:);DATA2(:)]);
    
    
    A = DATA;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential.DATA = A;
    
    A = DATA2;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential2.DATA = A;
    
    tmp = Differential.DATA;
    tmp = tmp(:,[4 8 9:12],:);
    tmp2 = Differential2.DATA;
    tmp2 = tmp2(:,[4 8 9:12],:);
    Max = max([tmp(:);tmp2(:)]);
    Min = min([tmp(:);tmp2(:)]);

    clear DATA DATA2 Idx tmp tmp2 A
    
    DATA = AllSubjects_Data(5).SubROI(iROI).DATA(:,:,Include);
    
    MEAN = nanmean(DATA,3);
    SEM = nansem(DATA,3);
    STD = nanstd(DATA,3);
    
    Name = [AllSubjects_Data(5).SubROI(iROI).name '_MNI ROI - With subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM);
    ToPlot.Subjects.Data = DATA;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
%     PlotLayers(ToPlot, Name, Visible, PRINT)

    clear ToPlot

    
    A = DATA;
    A = [A A(:,1:3,:)-A(:,4:6,:)];

    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    DATA = A;
    
    MEAN = nanmean(DATA,3);
    SEM = nansem(DATA,3);
    STD = nanstd(DATA,3);
    
    Name = [AllSubjects_Data(5).SubROI(iROI).name '_MNI ROI - Differential with subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM);
    ToPlot.Subjects.Data = DATA;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
    ToPlot.Max = Max;
    ToPlot.Min = Min;

    PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
end


for iROI = 1:7
    disp(AllSubjects_Data(6).SubROI(iROI).name)
    
    Idx = find(strcmp({AllSubjects_Data_Ori.name},...
        ['r' AllSubjects_Data(6).SubROI(iROI).name '_MNI.nii']));
    DATA = AllSubjects_Data_Ori(Idx).DATA(:,:,Include);
    DATA2 = AllSubjects_Data_Ori(Idx).DATA2(:,:,Include);
    MAX = max([DATA(:);DATA2(:)]);
    MIN = min([DATA(:);DATA2(:)]);
    
    A = DATA;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential.DATA = A;
    
    A = DATA2;
    if strcmp(ROIs{iROI}(1),'V')
        A = [A A(:,4:6,:)-A(:,1:3,:)];
    else
        A = [A A(:,1:3,:)-A(:,4:6,:)];
    end
    
    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    Differential2.DATA = A;
    
    tmp = Differential.DATA;
    tmp = tmp(:,[4 8 9:12],:);
    tmp2 = Differential2.DATA;
    tmp2 = tmp2(:,[4 8 9:12],:);
    Max = max([tmp(:);tmp2(:)]);
    Min = min([tmp(:);tmp2(:)]);

    clear DATA DATA2 Idx tmp tmp2 A

    
    DATA = AllSubjects_Data(6).SubROI(iROI).DATA(:,:,Include);
    
    MEAN = nanmean(DATA,3);
    SEM = nansem(DATA,3);
    STD = nanstd(DATA,3);
    
    Name = [AllSubjects_Data(6).SubROI(iROI).name '_MNI ROI - With subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM);
    ToPlot.Subjects.Data = DATA;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
%     PlotLayers(ToPlot, Name, Visible, PRINT)
    
    clear ToPlot
    
    
    A = DATA;

    A = [A A(:,4:6,:)-A(:,1:3,:)];

    for i=1:3
        TEMP(:,i,:) = A(:,i*3,:) - ( A(:,i*3-1,:)+A(:,i*3-2,:) );
    end
    A = [A(:,1:3,:) TEMP(:,1,:) A(:,4:6,:) TEMP(:,2,:) A(:,7:9,:) TEMP(:,3,:)];
    clear TEMP
    
    DATA = A;
    
    MEAN = nanmean(DATA,3);
    SEM = nansem(DATA,3);
    STD = nanstd(DATA,3);
    
    Name = [AllSubjects_Data(6).SubROI(iROI).name '_MNI ROI - Differential with subjects'];
    ToPlot.shadedErrorBar.mean = flipud(MEAN);
    ToPlot.shadedErrorBar.errorbar = flipud(SEM);
    ToPlot.Subjects.Data = DATA;
    ToPlot.MAX = MAX;
    ToPlot.MIN = MIN;
    ToPlot.Subjects.List = SubjectList(Include,:);
    ToPlot.Subjects.Legend.Bool = 0;
    ToPlot.Subjects.Legend.Legend = SubjectList(Include,:);
    ToPlot.Include=Include;
    
    ToPlot.Max = Max;
    ToPlot.Min = Min;

    PlotLayers(ToPlot, Name, Visible, PRINT)

    clear ToPlot    
end

for iROI=1:numel(ROIs)
    Name = [ROIs{iROI} ' SPM - Differential with subjects'];
    Command = [Command ' ' strrep(Name, ' ', '_') '.pdf']; %#ok<*AGROW>
    Name = [ROIs{iROI} ' ROI - Differential with subjects'];
    Command = [Command ' ' strrep(Name, ' ', '_') '.pdf'];
    Name = [ROIs{iROI} ' SyN - Differential with subjects'];
    Command = [Command ' ' strrep(Name, ' ', '_') '.pdf'];
end

system([...
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite ' ...
    '-sOutputFile=' fullfile(StartFolder, 'Coregistration_Comparisons.pdf') ' ' Command])