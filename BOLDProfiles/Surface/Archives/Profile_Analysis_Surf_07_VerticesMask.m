%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');

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

NbLayers = 6;
NbLayers = NbLayers+1;


for SubjInd = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    Data_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'BetaMapping','8Surf');
    
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Surfaces');
    
    % Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')

    Tab = [0;1;2];
    
    % For the 2 hemispheres
    for hs = 1:2
        
        if hs==1
            fprintf('   Left hemipshere\n')
            HsSufix = 'l';
        else
            fprintf('   Right hemipshere\n')
            HsSufix = 'r';
        end
        
        Betas = dir(fullfile(Data_Folder, ['Beta*' HsSufix 'cr.vtk']));
        
        [VertexBeta,FaceBeta,MappingBeta] = read_vtk(fullfile(Data_Folder,Betas(1).name),NbLayers, 1);
        NbVertices(hs,1)=size(VertexBeta,2);
        
        InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
            ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
        
        [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);
        NbVertices(hs,2)=size(Vertex,2);
        
        write_vtk(fullfile(Results_Folder,['Subj_' SubjID '_' HsSufix 'cr_ctrl_mask.vtk']), Vertex, Face, mean(MappingBeta))
        
        Mapping = zeros(size(Mapping));
        
        if NbVertices(hs,1)~=NbVertices(hs,2)
            NbVertices
            error('Number of vertices between beta mapping and source surface are different.')
        end

        Mapping = ones(size(Mapping))*2;
        Mapping(any(abs(MappingBeta)<.001,1))=1;
        Mapping(any(MappingBeta==0,1))=0;
        
        for iROI=1:numel(ROI)
            tmp = Mapping(ROI(iROI).VertOfInt{hs});
            ROI(iROI).Tab(:,1)=[0;.001;1]; %#ok<*SAGROW>
            ROI(iROI).Tab(1,hs+1)=sum(tmp==0)/numel(tmp);
            ROI(iROI).Tab(2,hs+1)=sum(tmp==1)/numel(tmp);
            ROI(iROI).Tab(3,hs+1)=sum(tmp==2)/numel(tmp);
        end
        
        tabulate(Mapping)
        T = tabulate(Mapping);
        Tab(:,hs+1) = T(:,3);
        
        write_vtk(fullfile(Results_Folder,['Subj_' SubjID '_' HsSufix 'cr_mask.vtk']), Vertex, Face, Mapping)
        
        clear Betas  Vertex Face Mapping VertexBeta FaceBeta MappingBeta
        
    end
    
    save(fullfile(Results_Folder,['Subj_' SubjID '_brain_coverage.mat']),'Tab')
    save(fullfile(SubjectFolder,'ROI_MIPAV',['Subj_' SubjID '_ROI_coverage.mat']),'ROI')

end

