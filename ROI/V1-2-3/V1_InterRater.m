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
%     '15';...
    '16'
    ];




for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    for hs=1:2
        %% Get inflated surface
        cd(fullfile(SubjectFolder, 'Structural', 'CBS'))
        if hs==1
            inf_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg_inf.vtk']);
            norm_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_lcr_gm_avg.vtk']);
        else
            inf_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_rcr_gm_avg_inf.vtk']);
            norm_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_rcr_gm_avg.vtk']);
        end
        
        [vertex,face,mapping] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', inf_vtk.name), 0, 1);
        [norm_vertex,norm_face,mapping_norm] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', norm_vtk.name), 0, 1);
        
        mapping = zeros(size(mapping));
        
        
        %% Get data
        cd(fullfile(StartFolder, 'Subjects_Data', 'ROI_V1'))
        if hs==1
            fprintf(' Left HS\n')
            suffix = 'l';
        else
            fprintf(' Right HS\n')
            suffix = 'r';
        end
        
        ROI_vtk = dir(['T1_' SubjID '_' suffix 'cr_V1.vtk']);
        [~,~,ROI1] = read_vtk(ROI_vtk.name, 0, 1);  
        unique(ROI1)
        VertOfInt1 = find(any([ROI1==22;ROI1==122]));
        mapping(VertOfInt1)= mapping(VertOfInt1) + 1;
        clear ROI_vtk
        
        ROI_vtk = dir(['T1_' SubjID '_' suffix 'cr_V1_un.vtk']);
        [~,~,ROI2] = read_vtk(ROI_vtk.name, 0, 1);    
        unique(ROI2)
        VertOfInt2 = find(any([ROI2==22;ROI2==122]));
        mapping(VertOfInt2)= mapping(VertOfInt2) + 2;
        
        DiceCoeff(SubjInd,hs)  = 2*numel(intersect(VertOfInt2,VertOfInt1)) / (numel(VertOfInt2)+numel(VertOfInt1));
        
%         M = [logical(ROI1)' logical(ROI2)'];
%         Vert = unique([VertOfInt2 VertOfInt1]);
%         M = M(Vert,:);
%         r = ICC(M, '1-1')

        if hs==1
            write_vtk(['Subj_' SubjID '_V1_lcr_RG_UN.vtk'], vertex, face, mapping)
            write_vtk(['Subj_' SubjID '_V1_lcr_norm_RG_UN.vtk'], norm_vertex, norm_face, mapping)
        else
            write_vtk(['Subj_' SubjID '_V1_rcr_RG_UN.vtk'], vertex, face, mapping)
            write_vtk(['Subj_' SubjID '_V1_rcr_norm_RG_UN.vtk'], norm_vertex, norm_face, mapping)
        end

        clear vertex face mapping ROI2 ROI1 Vert M
        
        
    end
    
    
end


cd(StartFolder)