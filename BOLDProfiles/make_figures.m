%% uses CSV files of saved data to plot the data for article

clear; close all; clc;

NbLayers = 6;
NbRuns = 4;
NbCdt = 6;
NbBlocks = 3;


DesMat = (1:NbLayers)-mean(1:NbLayers);
% DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

PlotSubjects = 0;
Switch = 1;
FigDim = [100 100 1800 1000];
Visible = 'on';
Transparent = 1;
Fontsize = 12;
Scatter = linspace(0,.4,11);


CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))
Get_dependencies('/home/remi')

FigureFolder = fullfile(CodeFolder, 'Figures', strcat(num2str(NbLayers+2), '_layers'));

% DataFolder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T';
DataFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

Results_Folder = fullfile(DataFolder, 'DataToExport');

suffix = {...
    '_stim-A_att-A';...
    '_stim-V_att-A';...
    '_stim-AV_att-A';...
    '_stim-A_att-V';...
    '_stim-V_att-V';...
    '_stim-AV_att-V';...
    };

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

NbSubj = size(SubjectList,1);

for iSubj=1:NbSubj
    sets{iSubj} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3'};

NbROI = size(ROIs,1);


%% Initialize variables.
filename = fullfile(Results_Folder, 'group_data-surf_ROI-A1_hs-both.csv');
delimiter = ',';
endRow = 783;

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%[^\n\r]';

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

clearvars filename delimiter endRow formatSpec fileID dataArray ans;

% cdt vectors (one column for each cdt)
CdtVec = false(size(RowName,1), NbCdt);
for iCdt = 1:NbCdt
    CdtVec(:, iCdt) = contains(RowName, suffix{iCdt});
end

SubjVec = false(size(RowName,1), NbSubj);
for iSubj = 1:NbSubj
    SubjVec(:, iSubj) = contains(RowName, ['sub-' SubjectList(iSubj,:)]);
end



V_stim = [2 5];
% Row2Select = any(CdtVec:V_stim);

for iSubj = 1:NbSubj
    
    Subj_Data = nan(NbRuns*NbBlocks, NbLayers, numel(V_stim));
    
    for iCdt = 1:numel(V_stim)
        Rows2Choose = [SubjVec(:,iSubj), CdtVec(:,V_stim(iCdt))];
        Rows2Choose = all(Rows2Choose, 2);
        Subj_Data(:,:,iCdt) = Data(Rows2Choose,:);
    end
    
    Subj_Data = mean(Subj_Data,3);
    
    All_Subjs_Profile(iSubj,:) = mean(Subj_Data,1);
    
    X=repmat(DesMat, [size(Subj_Data, 1), 1] );
    Subj_Data=Subj_Data';
    Subj_Data = Subj_Data(:);
    [B,~,~] = glmfit(X, Subj_Data, 'normal', 'constant', 'off');
    
    SubjectsBetas(iSubj,1:size(X,2)) = B;
end

% V_stim_Data = Data


DATA.WithSubj = PlotSubjects;

DATA.Scatter = Scatter;
DATA.WithPerm = 0;
DATA.PlotInset = 0;
DATA.YLabel = 'Param. est. [a u]';

close all

figure('position', FigDim, 'name', 'Deactivations', 'Color', [1 1 1], 'visible', Visible)

DATA.MVPA = 0;
DATA.YLabelInset = 1;
DATA.InsetLim = [.5 .4;-.5 -.4];
if DATA.WithSubj
    DATA.MIN = -0.6;
    DATA.MAX = 1;
else
    DATA.MIN = -0.3;
    DATA.MAX = 0.15;
end

DATA.OneSideTTest = {'left' 'both' 'both'};

Fontsize = 12;

subplot(2,1,1)
PlotRectangle(NbLayers,Fontsize,Switch)
subplot(2,1,1)


DATA.Name = char({'V vs. Baseline';'';'A1'});
DATA.Data = All_Subjs_Profile;
DATA.Betas = SubjectsBetas;
DATA.Color =  [0 0 0];
DATA.Thresholds = 0.05*ones(1,size(DATA.Betas,2));

PlotProfileAndBetas(DATA)

ax = subplot(2,1,2);
axis('off')
DATA.ax = ax.Position;
DATA.ToPermute = ToPermute;
PlotInsetFinal(DATA)


