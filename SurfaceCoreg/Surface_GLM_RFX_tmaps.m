%%
clear
clc

HS = 'L';
hs = 'l';
inf_hs = 'r';

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')));

SubjToInclude = true(13,1);
SubjToInclude([4 11],1) = false;

NbLayers = 6;

Conditions_Names = {...
    'A-Stim_Auditory-Attention', ...
    'V-Stim_Auditory-Attention', ...
    'AV-Stim_Auditory-Attention', ...
    'A-Stim_Visual-Attention', ...
    'V-Stim_Visual-Attention', ...
    'AV-Stim_Visual-Attention'};


for Smooth = 1
    
    DataFolder = fullfile('D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\SurfCoreg');
    
    InfSurfFile = fullfile(DataFolder, ['Surface_ls_' inf_hs 'h_inf.vtk']);
    [InfVertex,InfFace,InfMapping] = read_vtk(InfSurfFile, 0, 1);
    
    if Smooth
        suffix = '_smoothdata';
    else
        suffix = '';
    end
    
    
    %% Basic
    cd(fullfile(DataFolder, [HS 'H']))
    
    mkdir(fullfile(DataFolder, [HS 'H'], 'Baseline'))
    
    fprintf('\n Against baseline')
    for iCond = 1:6
        
        clear AllLayers AllLayersNoSmooth
        
        fprintf('\n Reading data\n')
        
        for iLayer = 1:NbLayers
            
            VTK_file = dir(['GrpSurf_' Conditions_Names{iCond} '_' num2str(iLayer) '-Layer_' hs 'h' suffix '.vtk']);
            
            disp(VTK_file.name)
            
            [Vertex,Face,Mapping] = read_vtk(VTK_file.name, 12, 1);
            
            AllLayers(:,:,iLayer) = Mapping; %#ok<SAGROW>
            
            if Smooth
                VTK_file = dir(['GrpSurf_' Conditions_Names{iCond} '_' num2str(iLayer) '-Layer_' hs 'h.vtk']);
                [~,~,Mapping] = read_vtk(VTK_file.name, 12, 1);
                AllLayersNoSmooth(:,:,iLayer) = Mapping; %#ok<SAGROW>
            end
            
        end
        
        AllLayers = AllLayers(SubjToInclude,:,:);
        
        if Smooth
            AllLayersNoSmooth = AllLayersNoSmooth(SubjToInclude,:,:);
            IsZero = AllLayersNoSmooth==0;
        else
            IsZero = AllLayers==0;
        end
        IsZero = any(IsZero,3);
        Vertex_to_include = sum(IsZero)<=(11-4);
        
        tmp = AllLayers;
        tmp(tmp==0)=nan;
        tmp =  nanmean(tmp,3);
        tmp = nanmean(tmp)./nansem(tmp);
        
        T_map = zeros(1, size(Mapping,2));
        T_map(Vertex_to_include) = tmp(Vertex_to_include);
        
        write_vtk(fullfile(DataFolder, [HS 'H'], 'Baseline', ...
            [Conditions_Names{iCond} '_' hs 'h_T_map_inf' suffix '.vtk']), InfVertex, InfFace, T_map')
        
        
    end
    
    
    %%
    clear AllLayers AllLayersNoSmooth
    
    for iCond = 1:6
        
        fprintf('\n Reading data\n')
        
        for iLayer = 1:NbLayers
            
            VTK_file = dir(['GrpSurf_' Conditions_Names{iCond} '_' num2str(iLayer) '-Layer_' hs 'h' suffix '.vtk']);
            
            disp(VTK_file.name)
            
            [Vertex,Face,Mapping] = read_vtk(VTK_file.name, 12, 1);
            
            AllLayers(:,:,iLayer,iCond) = Mapping;
            
        end
        
    end
    
    
    %% Average over attention
    fprintf('\n Average stim over attention')
    Avg = {[1;4], [2;5], [3;6]};
    Cdt = {'A_stim', 'V_stim', 'AV_stim'};
    
    for iCdt = 1:numel(Avg)
        
        AvgLayers = ...
            (sum(AllLayers(SubjToInclude,:,:,Avg{iCdt}(1,:)),4) ...
            + sum(AllLayers(SubjToInclude,:,:,Avg{iCdt}(2,:)),4))/2;
        
        tmp = AvgLayers;
        tmp(tmp==0)=nan;
        tmp =  nanmean(tmp,3);
        tmp = nanmean(tmp)./nansem(tmp);
        
        T_map = zeros(1, size(Mapping,2));
        T_map(Vertex_to_include) = tmp(Vertex_to_include);
        
        write_vtk(fullfile(DataFolder, [HS 'H'], 'Baseline', ...
            [Cdt{iCdt} '_' hs 'h_T_map_inf' suffix '.vtk']), InfVertex, InfFace, T_map')
        
        
    end
    
    
    %% Cross modal
    fprintf('\n Cross modal')
    mkdir(fullfile(DataFolder, [HS 'H'], 'CrossModal'))
    CrossModal = {...
        [3 6 ; 1 4], [3;1], [6;4],...
        [3 6 ; 2 5], [3;2], [6;5]};
    CrossModalName = {...
        'AV-A', 'AV-A_att_A', 'AV-A_att_V', ...
        'AV-V', 'AV-V_att_A', 'AV-V_att_V'};
    
    for iCrossModal = 1:numel(CrossModal)
        
        CrossModLayers = ...
            (sum(AllLayers(SubjToInclude,:,:,CrossModal{iCrossModal}(1,:)),4) ...
            - sum(AllLayers(SubjToInclude,:,:,CrossModal{iCrossModal}(2,:)),4))/numel(CrossModal{iCrossModal}(1,:));
        
        tmp = CrossModLayers;
        tmp(tmp==0)=nan;
        tmp =  nanmean(tmp,3);
        tmp = nanmean(tmp)./nansem(tmp);
        
        T_map = zeros(1, size(Mapping,2));
        T_map(Vertex_to_include) = tmp(Vertex_to_include);
        
        write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
            [CrossModalName{iCrossModal} '_' hs 'h_T_map_inf' suffix '.vtk']), InfVertex, InfFace, T_map')
        
    end
    
    
    %% Attention
    fprintf('\n Attention')
    mkdir(fullfile(DataFolder, [HS 'H'], 'Attention'))
    Att = {[1:3 ; 4:6], [1 ; 4], [2 ; 5], [3 ; 6]};
    AttName = {'Att_A-Att_V', 'Stim_A--Att_A-Att_V', 'Stim_V--Att_A-Att_V', 'Stim_AV--Att_A-Att_V'};
    
    for iAtt = 1:numel(Att)
        
        AttlLayers = ...
            (sum(AllLayers(SubjToInclude,:,:,Att{iAtt}(1,:)),4) ...
            - sum(AllLayers(SubjToInclude,:,:,Att{iAtt}(2,:)),4))/numel(CrossModal{iCrossModal}(1,:));
        
        tmp = AttlLayers;
        tmp(tmp==0)=nan;
        tmp =  nanmean(tmp,3);
        tmp = nanmean(tmp)./nansem(tmp);
        
        T_map = zeros(1, size(Mapping,2));
        T_map(Vertex_to_include) = tmp(Vertex_to_include);
        
        write_vtk(fullfile(DataFolder, [HS 'H'], 'Attention',...
            [AttName{iAtt} '_' hs 'h_T_map_inf' suffix '.vtk']), InfVertex, InfFace, T_map')
        
    end
    
    
end