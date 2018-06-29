%%
clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '16'
    ];




for SubjInd = size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    for hs=1:2
        %% Get inflated surface
        if hs==1
            fprintf(' Left HS\n')
            suffix = 'l';
        else
            fprintf(' Right HS\n')
            suffix = 'r';
        end
        
        cd(fullfile(SubjectFolder, 'Structural', 'CBS'))
        inf_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf.vtk']);
        norm_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg.vtk']);
        
        [vertex,face,mapping] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', inf_vtk.name), 0, 1);
        [norm_vertex,norm_face,mapping] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', norm_vtk.name), 0, 1);
        
        mapping = zeros(size(mapping));
        
        
        %% Get data
        cd(fullfile(StartFolder, 'Subjects_Data', 'ROI_TEs_cyt'))
        
        TEs_vtk = dir(['T1_' SubjID '_*' suffix 'cr_gm_avg_inf_TEs_Cyt.vtk']);
        
        [~,~,TEs] = read_vtk(TEs_vtk.name, 0, 1);
        
        
        cd(fullfile(StartFolder, 'Subjects_Data', 'ROI_TE'))
        
        TE_vtk = dir(['Subj_' SubjID '_' suffix 'cr_TE_RG.vtk']);
        [~,~,TE] = read_vtk(TE_vtk.name, 0, 1);
        
        
        %% Write
        VertOfInt1 = find(all([TEs==101;TE==1]));
        mapping(VertOfInt1)= 1;
        
        VertOfInt2 = find(all([TEs==102;TE==1]));
        mapping(VertOfInt2)= 2;
        
        VertOfInt3 = find(all([TEs==103;TE==1]));
        mapping(VertOfInt3)= 3;
        
        
        write_vtk(['Subj_' SubjID '_' suffix 'cr_inf_TEs.vtk'], vertex, face, mapping)
        write_vtk(['Subj_' SubjID '_' suffix 'cr_TEs.vtk'], norm_vertex, norm_face, mapping)
        
        
        clear vertex face mapping ROI2 ROI1 Vert M VertOfInt1 VertOfInt2 VertOfInt3
        
        
    end
    
    
end


cd(StartFolder)