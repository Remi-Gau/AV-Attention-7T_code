%%
% This script gets data from the beta values from the whole brain surface at the different depth
% for each condition, block, ROI, subject
clc; clear;

StartFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T';

Results_Folder = fullfile(StartFolder, 'DataToExport');
[~,~,~]=mkdir(Results_Folder);

addpath(genpath(fullfile(StartFolder, 'AV-Attention-7T_code', 'SubFun')))
Get_dependencies('D:\Dropbox\')

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


MinLayer = 1;
NbLayers = 7;
NbLayers = NbLayers+1;
Ind = NbLayers:-1:1;

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention'};


for SubjInd = 4:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    GLM_Folder = fullfile(SubjectFolder, 'FFX_Block');
    
    Data_Folder = fullfile(SubjectFolder,'BetaMapping','8Surf');
    
    % Creates a cell that lists the names of the beta images as well as
    % their column number in the design matrix
    load(fullfile(GLM_Folder, 'SPM.mat'))
    
    BetaNames = char(SPM.xX.name');
    
    % identify regressors of interests
    BetaOfInterest = ~any([BetaNames(:,9)=='T' BetaNames(:,7)=='R' BetaNames(:,7)=='c' ...
        strcmp(cellstr(BetaNames(:,end-1:end)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-2:end-1)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-3:end-2)), '2)') ...
        strcmp(cellstr(BetaNames(:,end-4:end-3)), '2)')], ...
        2);
    
    BetaOfInterest = find(BetaOfInterest);
    
    BetaNames = BetaNames(BetaOfInterest, :);
    
    % create vectors that identidy which condition and sessions of the data
    % matrix
    CdtVect = zeros(size(BetaNames,1), 1);
    for iCdt = 1:numel(Conditions_Names)
        row_idx = contains(cellstr(BetaNames), Conditions_Names{iCdt});
        CdtVect(row_idx) = iCdt;
    end
    
    SessVect = zeros(size(BetaNames,1), 1);
    for iSess = 1:4
        row_idx = contains( cellstr(BetaNames), sprintf('Sn(%i)', iSess) );
        SessVect(row_idx) = iSess;
    end
    
    clear SPM row_idx iSess iCdt
    
    
    %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID ...
        '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    
    %% Read features
    fprintf(' Reading data\n')
    
    % For the 2 hemispheres
    clear VertexWithDataHS MappingBothHS
    
    for hs = 1:2
        
        if hs==1
            fprintf('   Left hemipshere\n')
            HsSufix = 'l';
        else
            fprintf('   Right hemipshere\n')
            HsSufix = 'r';
        end
        
        FeatureSaveFile = fullfile(Data_Folder, ...
            [ 'Subj_' SubjID '_features_' HsSufix 'hs_' ...
            num2str(NbLayers) '_surf.mat']);
        
        %% Load data or extract them
        if exist(FeatureSaveFile, 'file')
            clear  AllMapping VertexWithData
            load(FeatureSaveFile, 'AllMapping', 'VertexWithData')
        else
            error('Data was not extracted from the VTK files.')
        end
        
        % initialize variables that stores data for the whole cortex but only
        % allocate data for vertices with have some data for
        Features = nan(NbVertex(hs),NbLayers,size(AllMapping,3));
        Features(VertexWithData,:,:) = AllMapping;
        
        fprintf(' Saving for ROI:\n')
        
        for iROI = 1:4
            
            fprintf(['  '  ROI(iROI).name '\n'])
            
            FileName = strcat(...
                'sub-', SubjID, ...
                '_data-surf_ROI-', ROI(iROI).name, ...
                '_hs-', HsSufix);
            
            Features_ROI = Features(ROI(iROI).VertOfInt{hs},:,:);
            
            % remove any vertices with one zero at any depth (CBS tools
            % gives a 0 value when no data was available)
            Vert2Rm = any( any(Features_ROI==0,3) ,2);
            Features_ROI(Vert2Rm,:,:) = [];
            
            Features_ROI = shiftdim(Features_ROI,2);
            
            % create a variable to know which column belongs to which depth
            LayerLabel = repmat(Ind, [size(Features_ROI, 2) 1]);
            LayerLabel = LayerLabel(:);
            
            Features_ROI = reshape(Features_ROI, ...
                [size(Features_ROI,1), size(Features_ROI, 2) * size(Features_ROI, 3)]);
            
            % saves to mat, csv and h5 format
            save(fullfile(Results_Folder, [FileName '.mat']), ...
                'Features_ROI', 'BetaNames', 'LayerLabel', ...
                'CdtVect', 'SessVect')
            
            csvwrite(fullfile(Results_Folder, [FileName '.csv']), ...
                Features_ROI)
            csvwrite(fullfile(Results_Folder, [FileName '_CdtVect' '.csv']), ...
                CdtVect)
            csvwrite(fullfile(Results_Folder, [FileName '_SessVect' '.csv']), ...
                SessVect)
            csvwrite(fullfile(Results_Folder, [FileName '_LayerLabel' '.csv']), ...
                LayerLabel)

            h5create(fullfile(Results_Folder, [FileName '.h5']),...
                ['/ROI-', ROI(iROI).name, '_hs-', HsSufix'], ...
                size(Features_ROI));
            h5write(fullfile(Results_Folder, [FileName '.h5']),...
                ['/ROI-', ROI(iROI).name, '_hs-', HsSufix'], ...
                Features_ROI)
            
            clear Features_ROI LayerLabel
            
        end
        
    end
    
    clear BetaNames CdtVect SessVect
    
end

