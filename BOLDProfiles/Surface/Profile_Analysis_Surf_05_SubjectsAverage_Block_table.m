%% saves each subject profile for each block and condition to a CSV file

clear; close all; clc;

NbLayers = 6;
NbRuns = 4;
NbCdt = 6;
NbBlocks = 3;
Median = 1;

StartFolder = fullfile(pwd, '..', '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder, 'Figures', 'ProfilesSurface', strcat(num2str(NbLayers), '_layers'));
[~,~,~] = mkdir(FigureFolder);

suffix = {...
    '_stim-A_att-A';...
    '_stim-V_att-A';...
    '_stim-AV_att-A';...
    '_stim-A_att-V';...
    '_stim-V_att-V';...
    '_stim-AV_att-V';...
    };
suffix = repmat(suffix', NbRuns*NbBlocks, 1);
suffix = suffix(:);
suffix =  char(suffix);
suffix = [repmat(num2str((1:3)'), NbRuns*NbCdt, 1) suffix];
suffix = [repmat('_block-', NbRuns*NbCdt*NbBlocks, 1) suffix];
tmp = repmat(1:4, NbBlocks, 1);
suffix = [repmat(num2str(tmp(:)), NbCdt, 1) suffix];
suffix = [repmat('_run-', NbRuns*NbCdt*NbBlocks, 1) suffix];


ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3'};

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

NbSubject = size(SubjectList,1);

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


%%
for iROI = 1:length(ROIs)
    AllSubjects_Data(iROI) = struct(...
        'name', ROIs{iROI}, ...
        'DATA', nan(NbLayers,6,NbSubject), ...
        'VertexCount', nan(NbSubject,NbLayers+1),...
        'Differential', struct('Blocks', struct('DATA', cell(1))));
end

%% Gets data for each subject
for SubjInd = 1:NbSubject
    
    cd(StartFolder)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf(['\n\n  Processing subject : ' SubjID '\n\n'])
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
        'Results', 'Profiles', 'Surfaces');
    
    
    for iROI=1:size(ROIs,1)
        
        File2Load = fullfile(SubjectFolder,...
            strcat('Data_Surf_Block_', ROIs{iROI}, '_', num2str(NbLayers+2), '_layers.mat'));
        
        if ~exist(File2Load, 'file')
            error('File %s is missing.', File2Load)
        else
            
            load(File2Load, 'Data_ROI')

            if isfield(Data_ROI, 'LayerMean')
                %% Compute main effect for each block
                if Median
                    Data = Data_ROI.LayerMedian(2:end-1,:,:);
                else
                    Data = Data_ROI.LayerMean(2:end-1,:,:); %#ok<*UNRCH>
                end
              Data = reshape(Data, [size(Data,2) * size(Data,3), size(Data,1)]);
              AllSubjects_Data(iROI).DATA{SubjInd} = Data;
                
            end
            
        end
        
        clear Data_ROIs SubROI_Ind Data_ROI
        
    end
    
end

clear SubjInd ROI_Ind

cd(StartFolder)