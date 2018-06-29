%%
clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
%     '03';...
%     '04';...
%     %     '06';...
%     '07';...
%     '08';...
%     '09';...
%     '11';...
%     '12';...
%     '13';...
%     %     '14';...
%     '15';...
%     '16'
    ];


for SubjInd = 1 %size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    cd(fullfile(StartFolder, 'Subjects_Data'))
    
    for hs=1
        
        if hs==1
            fprintf(' Doing left HS\n')
            Suffix = 'lcr';
        else
            fprintf(' Doing right HS\n')
            Suffix = 'rcr';
        end
        
        %% Get data
        cd(fullfile(StartFolder, 'Subjects_Data'))
        inf_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' Suffix '_gm_avg_inf.vtk']);
        [vertex,face,mapping] = read_vtk(inf_vtk.name, 0, 1);
        
        cd(fullfile(StartFolder, 'Subjects_Data'))
        qT1Map_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_qT1.vtk']);
        [~,~,qT1MapOri] = read_vtk(qT1Map_vtk.name, 0, 1);
        
        
        %% TE
        mapping = zeros(size(mapping));
        cd(fullfile(StartFolder, 'Subjects_Data'))
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te10.nii.vtk']);
        [~,~,TE10] = read_vtk(ROI_vtk.name, 0, 1);
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te11.nii.vtk']);
        [~,~,TE11] = read_vtk(ROI_vtk.name, 0, 1);
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te12.nii.vtk']);
        [~,~,TE12] = read_vtk(ROI_vtk.name, 0, 1);
         
        VertOfInt = find(any([TE10;TE11;TE12]>.1));
        T1Thresh = prctile(qT1MapOri(VertOfInt),50);
        mapping(VertOfInt) = 1;
        
        qT1Map = qT1MapOri;
        qT1Map = qT1Map<=T1Thresh;
        qT1Map = all([qT1Map; mapping]);
        VertOfInt = find(qT1Map);
        mapping(VertOfInt) = 2;
        
        write_vtk(['T1_' SubjID '_' Suffix  '_TE_p>0.1_T1_Thresh.vtk'], vertex, face, mapping)
        

        %% TE 1.0
%         mapping = zeros(size(mapping));
%         
%         VertOfInt = find(TE10>.1);
%         mapping(VertOfInt) = 1;
%         
%         qT1Map = qT1MapOri;
%         qT1Map = qT1Map<=T1Thresh;
%         qT1Map = all([qT1Map; mapping]);
%         VertOfInt = find(qT1Map);
%         mapping(VertOfInt) = 2;
%         
%         write_vtk(['T1_' SubjID '_' Suffix  '_TE1.0_p>0.1_T1_Thresh.vtk'], vertex, face, mapping)

        
        %% TE 1.1
%         mapping = zeros(size(mapping));
%         
%         VertOfInt = find(TE11>.1);
%         mapping(VertOfInt) = 1;
%         
%         qT1Map = qT1MapOri;
%         qT1Map = qT1Map<=T1Thresh;
%         qT1Map = all([qT1Map; mapping]);
%         VertOfInt = find(qT1Map);
%         mapping(VertOfInt) = 2;
%         
%         write_vtk(['T1_' SubjID '_' Suffix  '_TE1.1_p>0.1_T1_Thresh.vtk'], vertex, face, mapping)

        
        %% TE 1.2
%         mapping = zeros(size(mapping));
%         
%         VertOfInt = find(TE12>.1);
%         mapping(VertOfInt) = 1;
%         
%         qT1Map = qT1MapOri;
%         qT1Map = qT1Map<=T1Thresh;
%         qT1Map = all([qT1Map; mapping]);
%         VertOfInt = find(qT1Map);
%         mapping(VertOfInt) = 2;
%         
%         write_vtk(['T1_' SubjID '_' Suffix  '_TE1.2_p>0.1_T1_Thresh.vtk'], vertex, face, mapping)

        
         %% V1
%         mapping = zeros(size(mapping));
%         
%         cd(fullfile(StartFolder, 'Subjects_Data'))
%         ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmVisual_hOc1.nii.vtk']);
%         [~,~,V1] = read_vtk(ROI_vtk.name, 0, 1);
%         
%         VertOfInt = find(V1>.1);
%         mapping(VertOfInt) = 1;
%         
%         qT1Map = qT1MapOri;
%         qT1Map = qT1Map<=prctile(qT1Map(VertOfInt),50);
%         qT1Map = all([qT1Map; mapping]);
%         VertOfInt = find(qT1Map);
%         mapping(VertOfInt) = 2;
%         
%         write_vtk(['T1_' SubjID '_' Suffix  '_V1_p>0.1_T1_Thresh.vtk'], vertex, face, mapping)        

    end

    
end


cd(StartFolder)