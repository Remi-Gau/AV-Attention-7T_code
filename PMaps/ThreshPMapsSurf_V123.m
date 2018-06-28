%%
clc; clear; close all

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

Thresh = [10];

for SubjInd = 3 %size(SubjectList,1)
    
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
        [vertex,face,mapping] = read_vtk(inf_vtk.name, 0, 1);


        %% Get data
        cd(fullfile(SubjectFolder, 'PMap'))
 
        V1v_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi1.vtk']);        
        [~,~,V1v] = read_vtk(V1v_vtk.name, 0, 1);
        V1d_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi2.vtk']); 
        [~,~,V1d] = read_vtk(V1d_vtk.name, 0, 1);
        
        V1 = [V1v;V1d];
        
        V2v_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi3.vtk']);        
        [~,~,V2v] = read_vtk(V2v_vtk.name, 0, 1);
        V2d_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi4.vtk']); 
        [~,~,V2d] = read_vtk(V2d_vtk.name, 0, 1);
        
        V3v_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi5.vtk']);
        [~,~,V3v] = read_vtk(V3v_vtk.name, 0, 1);
        V3d_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_roi6.vtk']);
        [~,~,V3d] = read_vtk(V3d_vtk.name, 0, 1);
        
        V2_3 = [V2v;V2d;V3v;V3d];
        
        OtherROI = [];
        for iROI=[7:9 14:17]
            vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' ...
                suffix 'cr_gm_avg_inf_roi' num2str(iROI) '.vtk']);
            [~,~,tmp] = read_vtk(vtk.name, 0, 1);
            OtherROI(end+1,:) = tmp; %#ok<SAGROW>
        end

        SumV1 = sum(V1);
        SumV2_3 = sum(V2_3);
        SumOther = sum(OtherROI);

        for iThres=1:numel(Thresh)
            Mapping = zeros(2,size(V1v,2));
            
            Mapping(1,any(V1>Thresh(iThres)))=1;
            Mapping(1,SumV2_3>SumV1)=0;
            
            Mapping(2,any(V2_3>Thresh(iThres)))=2;
            Mapping(2,SumV2_3<SumV1)=0;
            Mapping(2,SumOther>SumV2_3)=0;
            
            Mapping = sum(Mapping);
            
            write_vtk(['Subj_' SubjID '_' suffix 'cr_V1_V2-3_thres' num2str(Thresh(iThres)) '.vtk'], vertex, face, Mapping)
        end
        
        clear vertex face mapping ROI2 ROI1 Vert M
        
        
    end
    
    
    %%
    FilesToZip = dir('rrw*.nii');
    for iFile=1:numel(FilesToZip)
        gzip(FilesToZip(iFile).name)
    end
    
    
end


cd(StartFolder)