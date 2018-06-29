clear
clc
close all


StartDirectory = pwd;

SubjectList = [...
    '02';...
%     '03';...
    '04';...
    %     '06';...
    '07';...
    '08';...
%     '09';...
    '11';...
    '12';...
    '13';...
    %     '14';...
    '15';...
    '16'
    ];

Visible = 'on';

FontSize = 6;


% Options for the SVM
opt.fs.do = 0; % feature selection
opt.rfe.do = 0; % recursive feature elimination
opt.permutation.test = 0;  % do permutation test
opt.session.curve = 0; % learning curves on a subsample of all the sessions

opt.print.do = 1;
opt.print.folder = fullfile(StartDirectory, 'Vertices_Analysis','Figures');

Analysis = struct('name', 'A Stim VS V Stim');
Analysis(end+1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');
Analysis(end+1) = struct('name', 'A Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'V Stim(A Att VS V Att)');
Analysis(end+1) = struct('name', 'AV Stim(A Att VS V Att)');


ROIs = {...
    'A1'; ... 1
    'TE1.0'; ...
    'TE1.1'; ...
    'TE1.2'; ...
    'PT'; ...
    
    'V1'; ... 6
    'V2'; ...
    'V3d'; ...
    'V3v'; ...
    'V4d'; ...
    'V4v'; ...
    'V5';...
    'V3';...
    'V4'};

Layers = [5 7];

%
SaveSufix = '_results_crosslayer';

if opt.fs.do
    SaveSufix = [SaveSufix '_FS'];
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

SaveSufix = [SaveSufix '.mat'];


Sides = {'LeftHS' ; 'RightHS' ; 'BothHS'};

NbLayers = [4 6];
NbLayers = NbLayers+1;

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\n\nAnalysing subject %s\n', SubjID)
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    SaveDir = fullfile(SubjectFolder, 'Analysis', 'Vertex', 'SVM');
    
    for iSVM=1:numel(Analysis)
        
        fprintf(' SVM: %s.\n', Analysis(iSVM).name)
        
        for iROI=1:size(ROIs,1)
            
            close all
            
            fprintf('  ROI: %s\n', ROIs{iROI,1})
            
            try
                
                load(fullfile(SaveDir, ...
                    ['SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} SaveSufix]), 'Results');
                
                Results = Results.session(end).rand(1).perm(1).CV;
                
                for iCV = 1:size(Results,1)
                    for iSide = 1:size(Results,2)
                        for iLayer=1:size(Results,3)
                            
                            grid = Results(iCV,iSide,iLayer).results{1}.gridsearch;
                            
                            if  numel(unique(grid.first.acc))>=2
                                name = ['GridSearch_Subj_' SubjID '_SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} ...
                                    '_CV_' Sides{iSide} '_All_Layers_' num2str(NbLayers(iLayer))]
                                GridSearchPlot(grid, name, opt)
                            end
                            
                            
                            
                            for iSubLayer=1:size(Results(iCV,iSide,iLayer).layers.results{1},2)
                                grid = Results(iCV,iSide,iLayer).layers.results{1}{iSubLayer}.grid;
                                
                                if  numel(unique(grid.first.acc))>=2
                                    name = ['GridSearch_Subj_' SubjID '_SVM_' Analysis(iSVM).name '_ROI_' ROIs{iROI} ...
                                        '_CV_' Sides{iSide} '_sublayer' num2str(iSubLayer) 'of' num2str(NbLayers(iLayer))]
                                    GridSearchPlot(grid, name, opt)
                                    
                                end
                                
                                
                            end
                        end
                    end
                end
                
                
            catch
            end
            clear Class_Acc iLayers
            
        end
        clear iROI
        
    end
    clear iSVM
    
end
clear NbRuns SubjInd