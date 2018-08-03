%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
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


MinLayer = 1;
NbLayers = 7;
NbLayers = NbLayers+1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};

NbWorkers = 4;

MatlabVer = version('-release');
if str2double(MatlabVer(1:4))>2013
    if isempty(gcp)
        KillGcpOnExit = 1;
        parpool(NbWorkers);
    else
        KillGcpOnExit = 0;
    end
    % else
    %     if matlabpool('size') == 0
    %         KillGcpOnExit = 1;
    %         matlabpool(NbWorkers)
    %     elseif matlabpool('size') ~= NbWorkers
    %         matlabpool close
    %         matlabpool(NbWorkers)
    %         KillGcpOnExit = 0;
    %     else
    %         KillGcpOnExit = 0;
    %     end
end

DesMat = (1:NbLayers-2)-mean(1:NbLayers-2);
% DesMat = [ones(NbLayers-2,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers-2,1) DesMat'];
DesMat = spm_orth(DesMat);


for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    GLM_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID], 'FFX_Block');
    
    Data_Folder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'BetaMapping','8Surf');
    
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Surfaces');
    [~,~,~]=mkdir(Results_Folder);
    
    mkdir(fullfile(Results_Folder,'Cdtions'))
    
    %% Creates a cell that lists the names of the beta images as well as their column number in the design matrix
    load(fullfile(GLM_Folder, 'SPM.mat'))
    
    BetaNames = char(SPM.xX.name');
    
    BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
        strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
        2);
    
    BetaOfInterest = find(BetaOfInterest);
    
    BetaNames(:,1:6)=[];
    
    clear SPM
    
    
        %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
    %% Read features
    fprintf(' Reading VTKs\n')
    
    % Format for reading the vertices from the VTK file
    Spec = repmat('%f ', 1, NbLayers);
    
    NbVertices = nan(1,2);
    
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
        
        FeatureSaveFile = fullfile(Data_Folder,[ 'Subj_' SubjID '_features_' HsSufix 'hs_' ...
            num2str(NbLayers) '_surf.mat']);
        
        InfSurfFile = fullfile(SubjectFolder, 'Structural','CBS', ...
            ['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' HsSufix 'cr_gm_avg_inf.vtk']);
        
        [Vertex,Face,Mapping] = read_vtk(InfSurfFile, 0, 1);
        NbVertices(hs)=size(Vertex,2);
        
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            load(FeatureSaveFile)
            VertexWithDataHS{hs} = VertexWithData; %#ok<*SAGROW>
        else
            
            AllMapping = nan(NbVertices(hs),NbLayers,size(Betas,1));
            
            fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
            
            parfor iBeta = 1:size(Betas,1)
                
                A = fileread(fullfile(Data_Folder,Betas(iBeta).name)); % reads file quickly
                B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
                
                C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
                Mapping = cell2mat(C); %clear C
                
                if size(Mapping,1)~=(NbVertices(hs)) %#ok<PFBNS>
                    error('A VTK file has wrong number of vertices:\n%s', fullfile(Data_Folder,Betas(iBeta).name))
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
        
        clear Betas InfSurfFile Mapping A B C
        

        %% Run GLMs for A, V and AV
        for CondInd = 1:2 % For each Condition
            
            Beta2Sel = [];
            Beta2Sel2 = [];
            for BlockInd = 1:3
                Beta2Sel = [Beta2Sel ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
                Beta2Sel2 = [Beta2Sel2 ;find(strcmp(cellstr(BetaNames), ...
                    [Conditions_Names{CondInd+3} ' - Block ' num2str(BlockInd) '*bf(1)']))]; %#ok<AGROW>
            end
            Beta2Sel = find(ismember(BetaOfInterest, Beta2Sel));
            Beta2Sel2 = find(ismember(BetaOfInterest, Beta2Sel2));
            
            Features = mean(cat(4,AllMapping(:,2:end-1,Beta2Sel), AllMapping(:,2:end-1,Beta2Sel2)),4);
            
            X=repmat(DesMat,size(Features,3),1);
            
            Y = shiftdim(Features,1);
            Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
            
            B = pinv(X)*Y;
            %             [B,~,~] = glmfit(X, Y, 'normal', 'constant', 'off');
            
            FullMapping = zeros(1,size(Vertex,2));
            FullMapping(VertexWithData) = B(1,:);
            
            
            
            Mapping = zeros(1,size(Vertex,2));
            
            nMax = 10;
            for iROI = 1:4
                Mapping = zeros(1,size(Vertex,2));
                tmp = FullMapping(ROI(iROI).VertOfInt{hs});
                tmp(tmp<0)=0;
                tmp(tmp>0)=1
                Mapping(ROI(iROI).VertOfInt{1}) = tmp;
                
                            write_vtk(fullfile(Results_Folder,'Cdtions',['Subj_' SubjID '_' HsSufix 'cr_' strrep(Conditions_Names{CondInd}(:,1:7),' ','') ...
                '_Cst_SortVert.vtk']), Vertex, Face, Mapping)
                
                
%                 Splits = floor(linspace(1,numel(tmp),nMax+1));
%                 [Y,I] = sort(tmp);
%                 
%                 for iSplit=1:nMax
%                     Mapping(ROI(iROI).VertOfInt{hs}(I(Splits(iSplit):Splits(iSplit+1))))=iSplit;
%                 end
%                 
% %                 Mapping(ROI(1).VertOfInt{1}) = tmp;
% 
%                 clear Y I tmp

            end
            
            
            write_vtk(fullfile(Results_Folder,'Cdtions',['Subj_' SubjID '_' HsSufix 'cr_' strrep(Conditions_Names{CondInd}(:,1:7),' ','') ...
                '_Cst_SortVert.vtk']), Vertex, Face, Mapping)
            

            
            clear Features Beta2Sel Beta2Sel2 B X Y Mapping
            
        end
        
        
        
        
    end
    
end

if str2double(MatlabVer(1:4))>2013
    if KillGcpOnExit
        delete(gcp)
    end
else
end

cd(StartFolder)