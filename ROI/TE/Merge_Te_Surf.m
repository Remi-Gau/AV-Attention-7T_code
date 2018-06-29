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
    '15';...
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
        
%         [vertex,face,mapping] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', inf_vtk.name), 0, 1);
        [norm_vertex,norm_face,mapping] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', norm_vtk.name), 0, 1);
        
        mapping = zeros(size(mapping));
        
        
        %% Get data
        cd(fullfile(StartFolder, 'Subjects_Data', 'ROI_TE'))
        
        
        ROI_vtk = dir(['T1_' SubjID '_' suffix 'cr_TE1.0_p>0.1_T1_Thresh_RG.vtk']);
        [~,~,TE10] = read_vtk(ROI_vtk.name, 0, 1);
        VertOfInt1 = find(any([TE10==35;TE10==135]));
        mapping(VertOfInt1)= mapping(VertOfInt1) + 1;
        
        ROI_vtk = dir(['T1_' SubjID '_' suffix 'cr_TE1.1_p>0.1_T1_Thresh_RG.vtk']);
        [~,~,TE11] = read_vtk(ROI_vtk.name, 0, 1);
        VertOfInt1 = find(any([TE11==35;TE11==135]));
        mapping(VertOfInt1)= mapping(VertOfInt1) + 2;
        
        ROI_vtk = dir(['T1_' SubjID '_' suffix 'cr_TE1.2_p>0.1_T1_Thresh_RG.vtk']);
        [~,~,TE11] = read_vtk(ROI_vtk.name, 0, 1);
        VertOfInt1 = find(any([TE11==35;TE11==135]));
        mapping(VertOfInt1)= mapping(VertOfInt1) + 4;
        
        
%         write_vtk(['Subj_' SubjID '_' suffix 'cr_inf_TE_p>0.1_T1_Thresh_RG.vtk'], vertex, face, mapping)
        write_vtk(['Subj_' SubjID '_' suffix 'cr_TE_RG.vtk'], norm_vertex, norm_face, mapping>1)
%         write_vtk(['Subj_' SubjID '_' suffix 'cr_TE_p>0.1_T1_Thresh_RG.vtk'], norm_vertex, norm_face, mapping)

        
        clear vertex face mapping ROI2 ROI1 Vert M
        
        
    end
    
    
end


cd(StartFolder)