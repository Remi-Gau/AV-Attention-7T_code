%% saves each subject profile for each block and condition to a CSV file

clear; close all; clc;

NbLayers = 6;
NbRuns = 4;
NbCdt = 6;
NbBlocks = 3;
Median = 1;

StartFolder = fullfile(pwd, '..', '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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

% Creates the following 
%     'run-1_block-1_stim-A_att-A '
%     'run-1_block-2_stim-A_att-A '
%     'run-1_block-3_stim-A_att-A '
%     'run-2_block-1_stim-A_att-A '
%     'run-2_block-2_stim-A_att-A '
%     'run-2_block-3_stim-A_att-A '
%     'run-3_block-1_stim-A_att-A '
%     'run-3_block-2_stim-A_att-A '
%     'run-3_block-3_stim-A_att-A '
%     'run-4_block-1_stim-A_att-A '
%     'run-4_block-2_stim-A_att-A '
%     'run-4_block-3_stim-A_att-A '
%     'run-1_block-1_stim-V_att-A '
%     'run-1_block-2_stim-V_att-A '
%     'run-1_block-3_stim-V_att-A '
%     'run-2_block-1_stim-V_att-A '
%     'run-2_block-2_stim-V_att-A '
%     'run-2_block-3_stim-V_att-A '
%     'run-3_block-1_stim-V_att-A '
%     'run-3_block-2_stim-V_att-A '
%     'run-3_block-3_stim-V_att-A '
%     'run-4_block-1_stim-V_att-A '
%     'run-4_block-2_stim-V_att-A '
%     'run-4_block-3_stim-V_att-A '
%     'run-1_block-1_stim-AV_att-A'
%     'run-1_block-2_stim-AV_att-A'
%     'run-1_block-3_stim-AV_att-A'
%     'run-2_block-1_stim-AV_att-A'
%     'run-2_block-2_stim-AV_att-A'
%     'run-2_block-3_stim-AV_att-A'
%     'run-3_block-1_stim-AV_att-A'
%     'run-3_block-2_stim-AV_att-A'
%     'run-3_block-3_stim-AV_att-A'
%     'run-4_block-1_stim-AV_att-A'
%     'run-4_block-2_stim-AV_att-A'
%     'run-4_block-3_stim-AV_att-A'
%     'run-1_block-1_stim-A_att-V '
%     'run-1_block-2_stim-A_att-V '
%     'run-1_block-3_stim-A_att-V '
%     'run-2_block-1_stim-A_att-V '
%     'run-2_block-2_stim-A_att-V '
%     'run-2_block-3_stim-A_att-V '
%     'run-3_block-1_stim-A_att-V '
%     'run-3_block-2_stim-A_att-V '
%     'run-3_block-3_stim-A_att-V '
%     'run-4_block-1_stim-A_att-V '
%     'run-4_block-2_stim-A_att-V '
%     'run-4_block-3_stim-A_att-V '
%     'run-1_block-1_stim-V_att-V '
%     'run-1_block-2_stim-V_att-V '
%     'run-1_block-3_stim-V_att-V '
%     'run-2_block-1_stim-V_att-V '
%     'run-2_block-2_stim-V_att-V '
%     'run-2_block-3_stim-V_att-V '
%     'run-3_block-1_stim-V_att-V '
%     'run-3_block-2_stim-V_att-V '
%     'run-3_block-3_stim-V_att-V '
%     'run-4_block-1_stim-V_att-V '
%     'run-4_block-2_stim-V_att-V '
%     'run-4_block-3_stim-V_att-V '
%     'run-1_block-1_stim-AV_att-V'
%     'run-1_block-2_stim-AV_att-V'
%     'run-1_block-3_stim-AV_att-V'
%     'run-2_block-1_stim-AV_att-V'
%     'run-2_block-2_stim-AV_att-V'
%     'run-2_block-3_stim-AV_att-V'
%     'run-3_block-1_stim-AV_att-V'
%     'run-3_block-2_stim-AV_att-V'
%     'run-3_block-3_stim-AV_att-V'
%     'run-4_block-1_stim-AV_att-V'
%     'run-4_block-2_stim-AV_att-V'
%     'run-4_block-3_stim-AV_att-V'


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


%%
for iROI = 1:length(ROIs)
    AllSubjects_Data(iROI) = struct(...
        'name', ROIs{iROI}, ...
        'DATA', []);
end

%% Gets data for each subject
for SubjInd = 1:NbSubject

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
                
                % Extract data for each layer / condition / block
                if Median
                    Data = Data_ROI.LayerMedian(2:end-1,:,:);
                else
                    Data = Data_ROI.LayerMean(2:end-1,:,:); %#ok<*UNRCH>
                end
                
                % Reshape to match the 'suffix' organization 
                Data = shiftdim(Data,1);
                Data = reshape(Data, [size(Data,1) * size(Data,2), size(Data,3)]);
                Data = fliplr(Data);
                
                AllSubjects_Data(iROI).DATA = cat(1, ...
                    AllSubjects_Data(iROI).DATA, ...
                    Data);

            end
            

            
        end
        
        clear Data_ROIs SubROI_Ind Data_ROI
        
    end
    
end

clear SubjInd ROI_Ind

prefix = repmat(cellstr(SubjectList)',[size(suffix,1), 1]);
prefix = [repmat('sub-', ) char(A(:)) ];


%                 Data = mat2cell(Data, ones(size(Data,1),1), ones(1,size(Data,2)));
%                 
%                 cat(2, cellstr(suffix), Data)
%             csvwrite(fullfile(Results_Folder, [FileName '.csv']), ...
%                 Features_ROI)
%             
%             fid = fopen (SavedTxt, 'w');
%             for iRow = 1:size(suffix,1)
%                 fprintf (fid, '%s,', Legends{i});
%             end

cd(StartFolder)