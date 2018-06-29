%%
clc; clear;

StartFolder=fullfile(pwd,'..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

ResultsFolder = fullfile(StartFolder,'Figures','8_layers');

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

Visibility='on';

ROIs(1) = struct('name', 'A1', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_V_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_V_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_A_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_A_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'A1_deact', 'NbVertices', []);

ROIs(end+1) = struct('name', 'PT', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_V_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_V_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_A_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_A_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'PT_deact', 'NbVertices', []);

ROIs(end+1) = struct('name', 'V1', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_V_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_V_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_A_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_A_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V1_deact', 'NbVertices', []);

ROIs(end+1) = struct('name', 'V2-3', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V2-3_V_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V2-3_V_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V2-3_A_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V2-3_A_deact', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V23_act', 'NbVertices', []);
ROIs(end+1) = struct('name', 'V23_deact', 'NbVertices', []);

COLOR_Subject = ColorSubject()

FigDim = [100 100 1800 1000];
Visible = 'off';



for SubjInd = 1:size(SubjectList,1)
    
    close all
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
   
    %% Load Vertices of interest for each ROI;
    load(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')

    for iROI = 1:numel(ROIs)
        Idx = find(strcmp(ROIs(iROI).name, {ROI.name}'));
        ROIs(iROI).NbVertices(SubjInd,1) = sum(cellfun('length', ROI(Idx).VertOfInt));
    end
    
end

for iROI = 1:numel(ROIs)
    ROIs(iROI).MinVert = floor(min(ROIs(iROI).NbVertices)/100)*100;
    MinVert = ROIs;
end
save(fullfile(ResultsFolder,'MinNbVert.mat'),'MinVert')


SubPlots = {[1 2 6 7], 3, 8, 4, 9, 5, 10};

for iROI=1:7:numel(ROI)

    
    %%
    figure('name', [ROI(iROI).name '- ROI'], 'Position', [100, 100, 1500, 1000], ...
        'Color', [1 1 1], 'Visible', Visibility);
    
    MAX = cat(2,ROIs(iROI:iROI+6).NbVertices);
    MAX = max(MAX(:));
    
    for i=iROI:iROI+6

        LEGEND={ROIs(i).name}';
        
        subplot(2,5,SubPlots{i-(iROI)+1})
        
        hold on
        
        errorbar(1,mean(ROIs(i).NbVertices), nansem(ROIs(i).NbVertices), '.k')
       
        for iSubj=1:size(SubjectList,1)
            plot(1+.1*iSubj,ROIs(i).NbVertices(iSubj), 'marker', '.', 'markersize', 30, ....
                'color', COLOR_Subject(iSubj,:))
        end
        
        axis([0.9 2.3 0 MAX])
        
        set(gca,'xtick',1+.1*(1:size(SubjectList,1)), 'xticklabel', SubjectList, 'ygrid', 'on',...
            'fontsize', 8)

        title(strrep(LEGEND, '_', '-'))

    end

    print(gcf, fullfile(ResultsFolder, [ROIs(iROI).name '_NbVertices.tif']), '-dtiff')
    
    
end

