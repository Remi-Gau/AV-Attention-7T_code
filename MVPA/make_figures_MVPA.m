%% uses CSV files of saved data to plot the data for the article

clear; close all; clc;

NbLayers = 6;

CodeFolder = 'D:\github\AV-Attention-7T_code';
DataFolder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T';
DependenciesFolder = 'D:';

% CodeFolder = '/home/remi/github/AV-Attention-7T_code';
% DataFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';
% DependenciesFolder = '/home/remi';

addpath(genpath(fullfile(CodeFolder, 'SubFun')))

FigureFolder = fullfile(CodeFolder, 'Figures', strcat(num2str(NbLayers+2), '_layers'));

Get_dependencies(DependenciesFolder)

Results_Folder = fullfile(DataFolder, 'DataToExport');

PlotSubjects = 1; % can be switched off (0) to no plot subject

% figure 3
% This will create some extra figures that are not in the paper (e.g AV-A
% in V1) but that can then act as positive controls
Cdt2Choose(1).name = '[AV VS A]_{att A, att V}';
Cdt2Choose(end).filename = 'Astim-vs-AVstim';
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'};
Cdt2Choose(2).name = '[AV VS V]_{att A, att V}';
Cdt2Choose(end).filename = 'Vstim-vs-AVstim';
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'};

% figure 4
Cdt2Choose(3).name = '[Att_A VS Att_V]_{A, V, AV}';
Cdt2Choose(end).test_side = {'both' 'both' 'both' 'both'};
Cdt2Choose(end).filename = 'AAtt-vs-VAtt';


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
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)']; % in case we want a quadratic component
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);


NbSubj = size(SubjectList,1);

% create permutations for exact sign permutation test
for iSubj=1:NbSubj
    sets{iSubj} = [-1 1]; %#ok<SAGROW>
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:}); clear sets
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];
clear a b c d e f g h i j k

NbROI = size(ROIs,1);

for iROI = 1:NbROI

    ROI_name = ROIs{iROI};

    for iCdt_2_plot = 1:numel(Cdt2Choose)

        %% get data
        filename = fullfile(Results_Folder, ...
            ['group_decoding_data-surf_ROI-' ...
            ROI_name '_Classification-' ...
            Cdt2Choose(iCdt_2_plot).filename...
            '_hs-both.csv']);

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

        % Subject vectors (one column for each subj)
        SubjVec = false(size(RowName,1), NbSubj);
        for iSubj = 1:NbSubj
            SubjVec(:, iSubj) = contains(RowName, ['sub-' SubjectList(iSubj,:)]);
        end
        clear iSubj


        %% does the math for each contrast
        for iSubj = 1:NbSubj

            % get data for that subject
            Rows2Choose = SubjVec(:, iSubj);
            Subj_Data = Data(Rows2Choose, :);

            % mean profile for that subject
            All_Subjs_Profile(iSubj, :) = mean(Subj_Data, 1); %#ok<SAGROW>

            % do laminar GLM
            X = repmat(DesMat, [size(Subj_Data, 1), 1] ); % design matrix

            % regorganize data
            Subj_Data=Subj_Data';
            Subj_Data = flipud(Subj_Data(:)-.5);

            [B,~,~] = glmfit(X, Subj_Data, 'normal', 'constant', 'off');

            SubjectsBetas(iSubj, 1:size(X, 2)) = B; %#ok<SAGROW>
        end


        % prepare for plotting
        DATA.WithSubj = PlotSubjects;
        DATA.FontSize = FontSize;
        DATA.Transparent = Transparent;
        DATA.YLabel = 'Param. est. [a u]';
        DATA.MVPA = 1;


        %% do actual plotting

        figure('position', FigDim, 'name', ' ', 'Color', [1 1 1], ...
            'visible', 'on')

        DATA.OneSideTTest = {Cdt2Choose(iCdt_2_plot).test_side{iROI} ...
            'both' 'both'};

        subplot(2, 1, 1)
        PlotRectangle(NbLayers, FontSize, Switch)
        subplot(2, 1, 1)


        DATA.Name = char({Cdt2Choose(iCdt_2_plot).name ; '' ; ROI_name});
        DATA.Data = fliplr(All_Subjs_Profile);
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
