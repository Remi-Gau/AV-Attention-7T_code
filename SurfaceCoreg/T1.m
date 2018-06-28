DataFolder = fullfile('/home/rxg243/Documents/GrpAverage/T1Avg');

InfSurfFile = fullfile(DataFolder, 'Surface_ls_lh_inf.vtk');
[InfVertex,InfFace,~] = read_vtk(InfSurfFile, 0, 1);

DataSurfFile = fullfile(DataFolder, 'Surface_ls_lh_trgsurf_groupavgdata.vtk');
[Vertex,Face,Mapping] = read_vtk(DataSurfFile, 0, 1);

write_vtk(fullfile(DataFolder, 'Grp_avg_T1_lh.vtk'), InfVertex, InfFace, Mapping')

InfSurfFile = fullfile(DataFolder, 'Surface_ls_rh_inf.vtk');
[InfVertex,InfFace,~] = read_vtk(InfSurfFile, 0, 1);

DataSurfFile = fullfile(DataFolder, 'Surface_ls_rh_trgsurf_groupavgdata.vtk');
[Vertex,Face,Mapping] = read_vtk(DataSurfFile, 0, 1);

write_vtk(fullfile(DataFolder, 'Grp_avg_T1_rh.vtk'), InfVertex, InfFace, Mapping')