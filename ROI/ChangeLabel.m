%%
clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SurfList = dir('T1_*cr_TE_SubDiv.vtk');

for iSurf = 2:numel(SurfList)
    
    disp(SurfList(iSurf).name)
    
    [vertex,face,mapping] = read_vtk(SurfList(iSurf).name, 0, 1);
    
    mapping(mapping==139)= 0;
    
%     mapping(mapping==102)= 20;
%     
%     mapping(mapping==103)= 30;
%     
%     c = find(all([mapping~=10;mapping~=20;mapping~=30]));
%     mapping(c)=0;
    
    write_vtk(SurfList(iSurf).name, vertex, face, mapping)
    
    clear vertex face mapping ROI2 ROI1 Vert M
    
end


cd(StartFolder)