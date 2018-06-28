%% Creates simple contrast at the groupe level
clear
clc

HS = 'L';
hs = 'l';
inf_hs = 'r';

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')));

Get_dependencies('/home/rxg243/Dropbox/')
Get_dependencies('D:\Dropbox')

NbLayers = 6;

Conditions_Names = {...
    'A-Stim_Auditory-Attention', ...
    'V-Stim_Auditory-Attention', ...
    'AV-Stim_Auditory-Attention', ...
    'A-Stim_Visual-Attention', ...
    'V-Stim_Visual-Attention', ...
    'AV-Stim_Visual-Attention'};

DataFolder = fullfile(StartFolder, 'SurfCoreg');

%% Stim

Stim = {[4 ; 1], [5 ; 2]};
StimName = {'Stim_A', 'Stim_V'};

for iStim = 1:numel(Stim)
    
    iStim
    Conditions_Names{Stim{iStim}(1)}
    VTK_file = fullfile(DataFolder, ['mean_inf_' Conditions_Names{Stim{iStim}(1)} '_lh.vtk']);
    [Vertex,Face,Mapping1] = read_vtk(VTK_file, 5, 1);
    Conditions_Names{Stim{iStim}(2)}
    VTK_file = fullfile(DataFolder, ['mean_inf_' Conditions_Names{Stim{iStim}(2)} '_lh.vtk']);
    [~,~,Mapping2] = read_vtk(VTK_file, 5, 1);

    Mapping = (Mapping1+Mapping2)/2;
    
    write_vtk(fullfile(DataFolder, ['mean_inf_' StimName{iStim} '_lh.vtk']),...
        Vertex, Face, Mapping(1:6,:)', 6)

end



%% Attention

Att = {[4 ; 1], [5 ; 2]};
AttName = {'Stim_A--Att_V-Att_A', 'Stim_V--Att_V-Att_A'};

for iAtt = 1:numel(Att)
    
    VTK_file = fullfile(DataFolder, ['mean_inf_' Conditions_Names{Att{iAtt}(1)} '_lh.vtk']);
    [Vertex,Face,Mapping1] = read_vtk(VTK_file, 5, 1);
    
    VTK_file = fullfile(DataFolder, ['mean_inf_' Conditions_Names{Att{iAtt}(2)} '_lh.vtk']);
    [~,~,Mapping2] = read_vtk(VTK_file, 5, 1);

    Mapping = Mapping1-Mapping2;
    
    write_vtk(fullfile(DataFolder, ['mean_inf_' AttName{iAtt} '_lh.vtk']),...
        Vertex, Face, Mapping', 6)

end

