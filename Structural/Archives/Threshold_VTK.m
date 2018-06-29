clear; clc;

addpath(fullfile(pwd,'SubFun'))

MainFolder = fullfile(pwd, 'Subjects_Data');

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

Smoothing = {'3' 'A';
             '6' 'B';
             '9' 'C';};

Thresholds = [1850 1900 1950];   
         
%  Folders definitions
for FWHM = 1:size(Smoothing,1)
    
    command = ['cp ' fullfile(MainFolder, 'Smooth', 'exp-00*', ['exp-00*-' Smoothing{FWHM,2}], ...
        'SmoothSurfaceMeshData', 'T1_*cr*.vtk ') fullfile(MainFolder)];
    system(command)
    
    for SubjInd = 5:size(SubjectList,1)
        
        % Subject's Identity
        SubjID = SubjectList(SubjInd,:);

        filename = dir(fullfile(MainFolder, ['T1_' SubjID '*lcr*data.vtk']));
        movefile(fullfile(MainFolder, filename.name), ...
            fullfile(MainFolder, ['T1_' SubjID '_lcr_' Smoothing{FWHM} '.vtk']))
        
        filename = dir(fullfile(MainFolder, ['T1_' SubjID '*rcr*data.vtk']));
        movefile(fullfile(MainFolder, filename.name), ...
            fullfile(MainFolder, ['T1_' SubjID '_rcr_' Smoothing{FWHM} '.vtk']))

        for iThresh = 1:numel(Thresholds)
            thres_vtk(fullfile(MainFolder, ['T1_' SubjID '_lcr_' Smoothing{FWHM} '.vtk']), 0, Thresholds(iThresh), [], 1, 1)
            thres_vtk(fullfile(MainFolder, ['T1_' SubjID '_rcr_' Smoothing{FWHM} '.vtk']), 0, Thresholds(iThresh), [], 1, 1)
        end
    end
    
end