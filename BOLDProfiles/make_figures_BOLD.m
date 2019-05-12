%% uses CSV files of saved data to plot the data for the article

clear; close all; clc;

NbLayers = 6;
NbRuns = 4;
NbCdt = 6;
NbBlocks = 3;

CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))

FigureFolder = fullfile(CodeFolder, 'Figures', strcat(num2str(NbLayers+2), '_layers'));

Get_dependencies('/home/remi')

% DataFolder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T';
DataFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

Results_Folder = fullfile(DataFolder, 'DataToExport');

PlotSubjects = 1; % can be switched off (0) to no plot subject

% figure 2
Cdt2Choose(1).name = 'A vs. Baseline';
Cdt2Choose(end).cdt = [1 4]; % auditory under A and V attention
Cdt2Choose(end).test_side = {'right' 'right' 'left' 'left'}; %side of permutation test
Cdt2Choose(2).name = 'V vs. Baseline';
Cdt2Choose(end).cdt = [2 5];
Cdt2Choose(end).test_side = {'left' 'left' 'right' 'right'};

% figure 3
% this will plot some extra contrast that are not in the paper (e.g AV-V
% for A1).
Cdt2Choose(3).name = '[AV - A]_{att A, att V}';
Cdt2Choose(end).cdt = [1 4 3 6]; % auditory under A and V attention ; AV under A and V attention ;
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'}; %side of permutation test
Cdt2Choose(4).name = '[AV - V]_{att A, att V}';
Cdt2Choose(end).cdt = [2 5 3 6]; % visual under A and V attention ; AV under A and V attention ;
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'};

% figure 4
% this will plot the results of A1 and PT upside down ([Att_V-Att_A] instead of ...
% [Att_A-Att_V]) compare to the figures in the paper
Cdt2Choose(5).name = '[Att_V - Att_A]_{A, V, AV}';
Cdt2Choose(end).cdt = [1 2 3 4 5 6]; % auditory under A and V attention ; AV under A and V attention ;
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'}; %side of permutation test


Transparent = 1;
Switch = 1;
FontSize = 12;
FigDim = [100 100 700 700];

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

suffix = {...
    '_stim-A_att-A';...
    '_stim-V_att-A';...
    '_stim-AV_att-A';...
    '_stim-A_att-V';...
    '_stim-V_att-V';...
    '_stim-AV_att-V';...
    };

% variable to read data file
delimiter = ',';
endRow = 800;

% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%[^\n\r]';


% Design matrix for laminar GLM
DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)']; % in case we want a
% quadratic component
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


NbSubj = size(SubjectList,1);

% create permutations for exact sign permutation test
for iSubj=1:NbSubj
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:}); clear sets
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
clear a b c d e f g h i j k

NbROI = size(ROIs,1);


for iROI = 1:NbROI
    
    ROI_name = ROIs{iROI};
    
    %% get data
    % Initialize variables.
    filename = fullfile(Results_Folder, ...
        ['group_data-surf_ROI-' ...
        ROI_name '_hs-both.csv']);

    % Open the text file.
    fileID = fopen(filename,'r');
    
    % Read columns of data according to the format.
    dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter, ...
        'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    % Close the text file.
    fclose(fileID);
    
    
    % Allocate imported array to column variable names
    RowName = dataArray{:, 1};
    Data = dataArray{:, 2};
    Data(:,2) = dataArray{:, 3};
    Data(:,3) = dataArray{:, 4};
    Data(:,4) = dataArray{:, 5};
    Data(:,5) = dataArray{:, 6};
    Data(:,6) = dataArray{:, 7};
    
    clearvars filename fileID dataArray
    
    % Condition vectors (one column for each cdt)
    CdtVec = false(size(RowName,1), NbCdt);
    for iCdt = 1:NbCdt
        CdtVec(:, iCdt) = contains(RowName, suffix{iCdt});
    end
    clear iCdt
    
    % Subject vectors (one column for each subj)
    SubjVec = false(size(RowName,1), NbSubj);
    for iSubj = 1:NbSubj
        SubjVec(:, iSubj) = contains(RowName, ['sub-' SubjectList(iSubj,:)]);
    end
    clear iSubj
    
    
    %% does the math for each contrast
    for iCdt_2_plot = 1:numel(Cdt2Choose)
        
        stim = Cdt2Choose(iCdt_2_plot).cdt;
        
        for iSubj = 1:NbSubj
            
            Subj_Data = nan(NbRuns*NbBlocks, NbLayers, numel(stim));
            
            % get data for that subject
            for iCdt = 1:numel(stim)
                Rows2Choose = [SubjVec(:, iSubj), CdtVec(:, stim(iCdt))];
                Rows2Choose = all(Rows2Choose, 2);
                Subj_Data(:, :, iCdt) = Data(Rows2Choose, :);
            end
            
            % in case we have to contrast between conditions 
            switch numel(stim)
                case 2
                case 4
                Subj_Data = Subj_Data(:, :, 3:4) - Subj_Data(:, :, 1:2);
                case 6
                Subj_Data = Subj_Data(:, :, 4:6) - Subj_Data(:, :, 1:3);
            end
            
            % mean across condition
            Subj_Data = mean(Subj_Data,3);
            
            % mean profile for that subject
            All_Subjs_Profile(iSubj, :) = mean(Subj_Data, 1); %#ok<SAGROW>
            
            % do laminar GLM
            X = repmat(DesMat, [size(Subj_Data, 1), 1] ); % design matrix
            
            % regorganize data
            Subj_Data=Subj_Data';
            Subj_Data = Subj_Data(:);
            
            [B,~,~] = glmfit(X, Subj_Data, 'normal', 'constant', 'off');
            
            SubjectsBetas(iSubj, 1:size(X, 2)) = B; %#ok<SAGROW>
        end
        
        
        % prepare for plotting
        DATA.WithSubj = PlotSubjects;
        DATA.FontSize = FontSize;
        DATA.Transparent = Transparent;
        DATA.YLabel = 'Param. est. [a u]';
        DATA.MVPA = 0;
        
        
        %% do actual plotting
        
        figure('position', FigDim, 'name', ' ', 'Color', [1 1 1], ...
            'visible', 'on')
        
        DATA.OneSideTTest = {Cdt2Choose(iCdt_2_plot).test_side{iROI} ...
            'both' 'both'};
        
        subplot(2, 1, 1)
        PlotRectangle(NbLayers, FontSize, Switch)
        subplot(2, 1, 1)
        
        
        DATA.Name = char({Cdt2Choose(iCdt_2_plot).name ; '' ; ROI_name});
        DATA.Data = All_Subjs_Profile;
        DATA.Betas = SubjectsBetas;
        DATA.Color =  [0 0 0];
        DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));
        
        PlotProfileAndBetas(DATA)
        
        ax = subplot(2, 1, 2);
        axis('off')
        DATA.ax = ax.Position;
        DATA.ToPermute = ToPermute;
        PlotInsetFinal(DATA)
        
    end
    
end