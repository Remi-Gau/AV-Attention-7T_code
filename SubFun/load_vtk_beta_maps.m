function [AllMapping, Face, Vertex, VertexWithData] = load_vtk_beta_maps(FeatureSaveFile, Data_Folder, HsSufix, NbVertices, NbLayers)
% loads data from VTK files or tries to extract them if the file is not
% there

%% Load data or extract them
if exist(FeatureSaveFile, 'file')
    load(FeatureSaveFile, 'AllMapping', 'Face', 'Vertex', 'VertexWithData')
else
    
    Betas = dir(fullfile(Data_Folder, ['Beta*' HsSufix 'cr.vtk']));
    
    AllMapping = nan(NbVertices(hs), NbLayers, size(Betas,1));
    
    fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
    
    parfor iBeta = 1:size(Betas,1)
        
        A = fileread(fullfile(Data_Folder, Betas(iBeta).name)); % reads file quickly
        B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
        
        C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
        Mapping = cell2mat(C); %clear C
        
        if size(Mapping,1)~=(NbVertices(hs)) %#ok<PFBNS>
            error('A VTK file has wrong number of vertices:\n%s', fullfile(Data_Folder, Betas(iBeta).name))
        end
        
        AllMapping(:,:,iBeta) = Mapping;
        
        fprintf(1,'\b.\n');
        
    end
    fprintf(1,'\b]\n');
    
    A = AllMapping==0;
    A = squeeze(any(A,2));
    A = ~any(A,2);
    VertexWithData = find(A);
    clear A
    
    AllMapping = AllMapping(VertexWithData,:,:);
    
    save(FeatureSaveFile,'Vertex','Face','AllMapping','VertexWithData', '-v7.3')
end

end