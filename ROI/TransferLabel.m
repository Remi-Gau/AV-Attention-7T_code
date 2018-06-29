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
    
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    for hs=1:2
        %% Get inflated surface
        if hs==1
            fprintf(' Left HS\n')
            suffix = 'l';
        else
            fprintf(' Right HS\n')
            suffix = 'r';
        end
    
    Surf = dir(['T1_' SubjID '_' suffix 'cr_TE_SubDiv.vtk']);
    [vertex,face,mapping] = read_vtk(Surf.name, 0, 1);
    
    cd ..
    InfSurf = dir(['T1_' SubjID '*' suffix 'cr_gm_avg_inf.vtk']);
    [infvertex,infface,~] = read_vtk(InfSurf.name, 0, 1);
    unique(mapping)
    
    new_mapping = zeros(size(mapping));
    new_mapping(mapping==4)= 4;
    new_mapping(mapping==5)= 5;
    new_mapping(mapping==6)= 6;
    
    cd('ROI_TEs_cyt')
    write_vtk(['T1_' SubjID '_' suffix 'cr_inf_TE_SubDiv.vtk'], infvertex, infface, mapping)
    
    clear vertex face mapping ROI2 ROI1 Vert M
    
    end
    
end


cd(StartFolder)