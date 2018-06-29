%%
clc; clear;

StartFolder=fullfile(pwd,'..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Structural','T1_profiles');

NbLayers = 11;

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

for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    Data_Folder = fullfile(SubjectFolder, 'Structural', 'CBS', 'Vertex', 'T1');
    
    Results_Folder = fullfile(SubjectFolder, 'Structural', 'CBS');
    
    
    %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
    %% Read features
    fprintf(' Reading T1 mappings\n')
    
    NbVertices = nan(1,2);
    
    % For the 2 hemispheres
    for hs = 1:2
        
        if hs==1
            fprintf('   Left hemipshere\n')
            HsSufix = 'l';
        else
            fprintf('   Right hemipshere\n')
            HsSufix = 'r';
        end
        
        T1MapSurfFile = dir(fullfile(Data_Folder, ...
            ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_' num2str(NbLayers) '_layers.vtk']));
        
        if exist(fullfile(Data_Folder,T1MapSurfFile.name), 'file')
            [Vertex,Face,Mapping] = read_vtk(fullfile(Data_Folder,T1MapSurfFile.name), NbLayers, 1);
            NbVertices(hs)=size(Vertex,2);
            AllMapping{hs} = Mapping;
        else
            error('File %s does not exist.',T1MapSurfFile)
        end
        
        clear Mapping
        
    end
    
    if any(NbVertex ~= NbVertices)
        NbVertex
        NbVertices
        error('The number of vertices does not match.')
    end
    
    for iROI=1:numel(ROI)
        
        disp(ROI(iROI).name)
        
        tmp  = [AllMapping{1}(1:end,ROI(iROI).VertOfInt{1}) ...
            AllMapping{2}(1:end,ROI(iROI).VertOfInt{2})];
        
        qT1_Profile(:,1) = nanmean(tmp,2);
        qT1_Profile(:,2) = nanstd(tmp,2);
        
        save(fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(iROI).name, '_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
        
        clear qT1_Profile
        
        %Left HS
        fprintf('   Left hemipshere\n')
        tmp  = AllMapping{1}(1:end,ROI(iROI).VertOfInt{1});
        qT1_Profile(:,1) = nanmean(tmp,2)
        qT1_Profile(:,2) = nanstd(tmp,2);
        save(fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(iROI).name, '_lh_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
        
        clear qT1_Profile
        
        %Right HS
        fprintf('   Right hemipshere\n')
        tmp  = AllMapping{2}(1:end,ROI(iROI).VertOfInt{2});
        qT1_Profile(:,1) = nanmean(tmp,2)
        qT1_Profile(:,2) = nanstd(tmp,2);
        save(fullfile(SubjectFolder,'Structural','CBS',...
            strcat('qT1_profile_', ROI(iROI).name, '_rh_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
        
        
        clear tmp qT1_Profile
    end
    
    
    tmp  = AllMapping{1}(1:end,:);
    qT1_Profile(:,1) = nanmean(tmp,2);
    qT1_Profile(:,2) = nanstd(tmp,2);
    clear tmp
    save(fullfile(SubjectFolder,'Structural','CBS',...
        strcat('qT1_profile_Cortex_lh_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
    
    
    tmp  = AllMapping{2}(1:end,:);
    qT1_Profile(:,1) = nanmean(tmp,2);
    qT1_Profile(:,2) = nanstd(tmp,2);
    clear tmp
    save(fullfile(SubjectFolder,'Structural','CBS',...
        strcat('qT1_profile_Cortex_rh_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
    
    
    tmp  = [AllMapping{1}(1:end,:) ...
        AllMapping{2}(1:end,:)];
    qT1_Profile(:,1) = nanmean(tmp,2);
    qT1_Profile(:,2) = nanstd(tmp,2);
    clear tmp
    save(fullfile(SubjectFolder,'Structural','CBS',...
        strcat('qT1_profile_Cortex_surf_', num2str(NbLayers), '.mat')), 'qT1_Profile')
    
    clear AllMapping qT1_Profile
    
    
end

cd(StartFolder)