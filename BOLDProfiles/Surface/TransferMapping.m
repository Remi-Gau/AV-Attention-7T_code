% scrit to transfer mapping from a vtk surface to its inflated counterpart

clc
clear
close all

StartDirectory = fullfile(pwd, '..', '..', '..');

addpath(genpath(fullfile(StartDirectory, 'SubFun')))
Get_dependencies('D:\Dropbox')

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

Source_dir = 'C:\Users\Remi\Documents';

for  iSubj = 1:size(SubjectList,1)
    SubjID = SubjectList(iSubj,:);
    
    T1 = spm_select('FPList', Source_dir, ...
        ['^T1_' SubjID '.*_lcr_gm_avg_inf_qT1.vtk$']);
    
    Files2transfer = spm_select('FPList', Source_dir, ...
        ['^Subj_' SubjID '_lcr_.*.vtk$']);

    
    [inf_vertex,inf_faces,T1] = read_vtk(T1, 0, 1);
    
    for  iFile = 1:size(Files2transfer,1)
        
        Filename = deblank(Files2transfer(iFile,:));
        
        [inf_vertex_func,inf_faces_func,mapping] = read_vtk(Filename, 0, 1);
        
        write_vtk([Filename(1:end-4) '_remap.vtk'], inf_vertex, inf_faces, mapping)
    end
end