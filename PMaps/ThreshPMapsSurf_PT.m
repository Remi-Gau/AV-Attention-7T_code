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

Surf2Load = {...
'IPL_6_6';...
'PoG_4_2';...
'PrG_6_5';...
'STG_6_1';...
'STG_6_2';...
'STG_6_3';...
'STG_6_4';...
'STG_6_5';...
'STG_6_6'};

Thresh = [40];

for SubjInd = 11%1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    for hs=2
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
 
        PT_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' suffix 'cr_gm_avg_inf_' Surf2Load{5} '.vtk']); 
        PT_vtk.name
        [~,~,PT] = read_vtk(PT_vtk.name, 0, 1);
        
%         OtherROI = [];
%         for iROI=[1:4 6:9]
%             vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' ...
%                 suffix 'cr_gm_avg_inf_' Surf2Load{iROI} '.vtk']);
%             [~,~,tmp] = read_vtk(vtk.name, 0, 1);
%             OtherROI(end+1,:) = tmp; %#ok<SAGROW>
%         end
%         SumOther = sum(OtherROI);

        [~,~,A1] = read_vtk(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1', ['Subj_' SubjID '_A1_' suffix 'cr_norm_RG_UN.vtk']), 0, 1);    

        for iThres=1:numel(Thresh)
            Mapping = zeros(1,size(PT,2));
            
            Mapping(1,PT>Thresh(iThres))=1;
%             Mapping(1,SumOther>PT)=0;
            Mapping(1,A1>0)=0;
            
            write_vtk(['Subj_' SubjID '_' suffix 'cr_PT_thres' num2str(Thresh(iThres)) '.vtk'], vertex, face, Mapping)
        end
        
        clear vertex face mapping ROI2 ROI1 Vert M
        
        
    end
    
    
    %%
%     FilesToZip = dir('rrw*.nii');
%     for iFile=1:numel(FilesToZip)
%         gzip(FilesToZip(iFile).name)
%     end
    
    
end


cd(StartFolder)