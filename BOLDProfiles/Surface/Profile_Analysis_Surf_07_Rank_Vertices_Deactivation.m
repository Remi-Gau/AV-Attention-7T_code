%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
%     '02';...
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

for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);

    
    %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    

    %% Read data
    fprintf(' Reading VTKs\n')

    NbVertices = nan(1,2);
    
    % For the 2 hemispheres
    for hs = 1
        
        if hs==1
            fprintf('   Left hemipshere\n')
            HsSufix = 'l';
        else
            fprintf('   Right hemipshere\n')
            HsSufix = 'r';
        end

        %%
        InfSurfFile = fullfile('C:\Users\Remi\Documents', ...
            ['Subj_' SubjID '_' HsSufix 'cr_AStim_Cst_smoothdata_remap.vtk']);
        
        [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);

        V1 = Mapping(ROI(3).VertOfInt{hs});
        [~,I] = sort(V1);
        I = I/numel(I);
        tmp = zeros(size(I));
        tmp(I<=.25) = -1;
        tmp(I>=.75) = 1;
        
        Mapping = zeros(1,size(Vertex,2));
        Mapping(ROI(3).VertOfInt{hs}) = tmp;
        write_vtk(fullfile('C:\Users\Remi\Documents',['Subj_' SubjID '_' HsSufix 'cr_AStim_Cst_sorted.vtk']), ...
            Vertex, Face, Mapping)
        
        clear Features Beta2Sel Beta2Sel2 B X Y Mapping

    end
    
end

cd(StartFolder)