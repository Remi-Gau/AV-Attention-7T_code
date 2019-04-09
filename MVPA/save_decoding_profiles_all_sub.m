clc; clear; close all;

NbLayersMVPA = 6;

% DataFolder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T';
DataFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

Results_Folder = fullfile(DataFolder, 'DataToExport');

% addpath(genpath(fullfile(CodeFolder, 'SubFun')))
% Get_dependencies('D:')

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

ROIs = {...
    'A1';...
    'PT';...
    'V1';...
    'V2-3'};

Analysis(1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');

opt.acroslayer.do = 0;
opt.leave2out.do = 0;

opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

opt.scaling.img.eucledian = 0;
opt.scaling.img.zscore = 0;
opt.scaling.feat.mean = 1;
opt.scaling.feat.range = 0;
opt.scaling.feat.sessmean = 0;
opt.scaling.idpdt = 1;

FFX = {'0'};

%% Set things
% init var
for iROI = 1:length(ROIs)
    for iSVM=1:numel(Analysis)
        AllSubjects_Data(iROI,iSVM) = struct(...
            'name', ROIs{iROI}, ...
            'analysis', Analysis(iSVM).name, ...
            'DATA', []);
    end
end

SaveSufix = '_results_surf_FixedC';
if opt.fs.do
    SaveSufix = [SaveSufix '_FS']; %#ok<*AGROW>
end
if opt.rfe.do
    SaveSufix = [SaveSufix '_RFE'];
end
if opt.permutation.test
    SaveSufix = [SaveSufix '_Perm'];
end
if opt.session.curve
    SaveSufix = [SaveSufix '_Lear'];
end
if opt.scaling.idpdt
    SaveSufix = [SaveSufix '_Idpdt'];
end
if opt.scaling.img.zscore
    SaveSufix = [SaveSufix '_ZScore'];
end
if opt.scaling.img.eucledian
    SaveSufix = [SaveSufix '_Eucl'];
end
if opt.scaling.feat.mean
    SaveSufix = [SaveSufix '_MeanCent'];
end
if opt.scaling.feat.range
    SaveSufix = [SaveSufix '_Range'];
end
if opt.scaling.feat.sessmean
    SaveSufix = [SaveSufix '_SessMeanCent'];
end

SaveSufix = [SaveSufix '_FWHM_' FFX{1} '_Layers_' num2str(NbLayersMVPA+2) '.mat'];


%% Get Data for MVPA
for SubjInd = 1:NbSubject
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf([' Processing subject : ' SubjID '\n'])
    
    %%
    for iSVM=1:numel(Analysis)
        
        for iROI=1:numel(ROIs)
            
            Save_vol = [...
                'SVM_' Analysis(iSVM).name...
                '_ROI_' ROIs{iROI} SaveSufix];
            
            load(fullfile(DataFolder, 'Subjects_Data', ['Subject_' SubjID],...
                'Transfer', 'SVM', Save_vol));
            
            CV = Results.session(end).rand.perm.CV;
            
            for iLayer= 2:NbLayersMVPA+1
                Label = [];
                Pred = [];
                for iCV=1:numel(CV)
                    Label(:,end+1) = CV(iCV).layers.results{1}{iLayer}.label;
                    Pred(:,end+1) = CV(iCV).layers.results{1}{iLayer}.pred(:,iLayer);
                end
                Acc(iLayer-1,:) = mean(Pred==Label,2)';
            end
            
            AllSubjects_Data(iROI,iSVM).DATA = cat(1, ...
                AllSubjects_Data(iROI,iSVM).DATA, ...
                Acc');
            
            clear Acc Pred Label
        end
        
    end
    
end

%% saves the data

NbRuns = 4;

for iSVM=1:size(AllSubjects_Data, 2)
    
    for iROI = 1:size(AllSubjects_Data, 1)
        
        Analysis = strrep(AllSubjects_Data(iROI,iSVM).analysis, 'VS', '-vs-');
        Analysis = strrep(Analysis, 'Stim', 'stim');
        Analysis = strrep(Analysis, ' ', '');
        
        FileName = ['group_decoding_data-surf_ROI-', ...
            AllSubjects_Data(iROI,iSVM).name, ...
            '_Classification-', Analysis,...
            '_hs-both'];
        
        Data = AllSubjects_Data(iROI,iSVM).DATA;
        
        
        % creates a label for each row
        NbBlocks = size(Data,1) / (NbSubject);
        suffix = repmat('_CV', [size(Data,1), 1]);
%         tmp = repmat(1:NbRuns, [NbBlocks,1]);
%         suffix = [suffix num2str(tmp(:))];
%         suffix = [suffix repmat('_block-', [NbRuns*NbBlocks,1])];
%         tmp = repmat((1:NbBlocks)', [NbRuns,1]);
%         suffix = [suffix num2str(tmp)];
        clear tmp
        
        prefix = repmat(cellstr(SubjectList)',[NbBlocks, 1]);
        prefix = char(prefix(:));
        prefix = [repmat('sub-', [size(prefix,1), 1] ) prefix ];
        
        labels = [prefix suffix];
        clear prefix suffix
        

        % save to .mat
        save(fullfile(Results_Folder, [FileName '.mat']), ...
            'Data', 'labels')
        
        
        % save to .csv
        fid = fopen (fullfile(Results_Folder, [FileName '.csv']), 'w');
        for iRow = 1:size(labels,1)
            fprintf (fid, '%s,', labels(iRow,:));
            fprintf (fid, '%f,', Data(iRow,:));
            fprintf (fid, '\n');
        end
        fclose (fid);
        
        clear Data
    end
    
end

