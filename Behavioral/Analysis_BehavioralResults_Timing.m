clear; clc; close all

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

FoldersNames = {...
    1:4;...
    1:4;...
    1:4;...
    1:2;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:4;...
    1:2;...
    1:4;...
    1:4;...
    };

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention', ...
    'A Targets - Auditory Attention', ...
    'A Targets - Visual Attention', ...
    'V Targets - Auditory Attention', ...
    'V Targets - Visual Attention', ...
    'Responses - Auditory Attention', ...
    'Responses - Visual Attention',...
    };

% Color for Conditions
COLOR =   [...
    255 0 0; ...
    255 100 100; ...
    255 200 200; ...
    0 255 0; ...
    100 255 100; ...
    200 255 200];
COLOR=COLOR/255;


IndStart = 5; % first row of data points in txt file

TargetTimeOut = 5; % s

%  Folders definitions
% StartFolder = fullfile(pwd,'..','..');
StartFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2');
addpath(genpath(fullfile(pwd,'..','..', 'SubFun')))


%% Define folders, number of runs, scans per run...
for SubjInd = 1:size(SubjectList,1)

    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    
%     figure('Name', ['Subject ' SubjID], 'Position', [100, 100, 1500, 1000])
    
    fprintf('\nRunning subject %s.\n\n', SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    BehavioralFolder = fullfile(SubjectFolder, 'Behavioral', 'Runs');
        
    NbRuns = length(FoldersNames{SubjInd});
    
    NbVol = nan(1,NbRuns);
    
    %% Collects the SOTs
    SOT_Blocks = cell(3,2,size(NbRuns,1));
    SOT_Targets = cell(2,2,size(NbRuns,1));
    SOT_Responses = cell(1,size(NbRuns,1));
    
    All_SOTs = cell(length(Conditions_Names),NbRuns);
    
    All_Durations = All_SOTs;
    
    fprintf('\nCollects the SOTs.\n\n')
    
    cd(BehavioralFolder)
    
    Block = ones(3,2);
    RespInBlock = nan(3,2,3);
    
    for RunInd = 1:NbRuns
        
        Block_Start = cell(3,2);
        
        Block_Duration = cell(3,2);
        
        Target_Index = cell(2,2);
        
        Response_Index = cell(2,2);
        
        Extra_Response = cell(1,1);

        %% Identify the relevant file and opens it
        tmp = [];
        EOF = [] ;
        fid = [];
        IndexStartExperiment =[];
        FileContent = {};
        TEMP = [];
        
        LogFileList = dir(strcat(...
            'Logfile_Subject_',  num2str(str2double(SubjID)), ...
            '_Run_', num2str(FoldersNames{SubjInd}(RunInd)), '*.txt'));
        disp(LogFileList.name)
        fid = fopen(fullfile (pwd, LogFileList.name));
        FileContent = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t', 'headerlines', IndStart, 'returnOnError',0);
        fclose(fid);
        
        
        % Finds the end of the file
        EOF = find(strcmp('Final_Fixation', FileContent{1,3}));
        if isempty(EOF)
            EOF = find(strcmp('Quit', FileContent{1,1}));
        end
        
        % Identify the time of the start of the XP
        IndexStartExperiment = find(strcmp('Pulse', FileContent{1,2}));
        IndexStartExperiment = sort(IndexStartExperiment);
        IndexStartExperiment = IndexStartExperiment(1);
        
        StartTime = str2double(FileContent{1,4}(IndexStartExperiment));
        
        % Counts the number of volumes and removed triggers from list
        TEMP = find(strcmp('30', FileContent{1,3}));
        NbVol(RunInd) = length(TEMP);
        
        Stim_Time = {FileContent{1,3}(IndexStartExperiment:EOF)  ...
            char(FileContent{1,4}(IndexStartExperiment:EOF))};
        
        TEMP = [ ...
            find(strcmp('SOA_Fix', Stim_Time{1,1})) ;...
            find(strcmp('SOA_Var', Stim_Time{1,1})) ;...
            find(strcmp('30', Stim_Time{1,1})) ;...
            find(strcmp('Fixation_Onset_Fix', Stim_Time{1,1})) ;...
            find(strcmp('Long_Fixation', Stim_Time{1,1}))];
        
        Stim_Time{1,1}(TEMP,:) = [];
        Stim_Time{1,2}(TEMP,:) = [];
        
        
        for i=1:size(Stim_Time{1,2},1)
            tmp(i,1)=str2double(Stim_Time{1,2}(i,:)); %#ok<SAGROW>
        end
        
        Stim_Time{1,2} = (tmp-StartTime)/10000;
        
        
        %% Identify targets and responses
        CurrentCondition = [];
        InBlock = 0;
        BlockType = [];
        BlockEnd = [];  %#ok<*NASGU>
        
        StimCount = 0;
        
        StimPresented = 0;
        TargetType = 0;
        Target = [];
        TargetPresentationTime = 0;
        
        BlockOrder = [];
        
        for i = 1:size(Stim_Time{1,2},1)
            
            if StimPresented && Stim_Time{1,2}(i)-TargetPresentationTime>TargetTimeOut
                StimPresented = 0;
                if Target==0
                    Target_Index{TargetType,CurrentCondition}(end+1)=TargetPresentationTime;
                else
                    if CurrentCondition==1
                        Target_Index{TargetType,2}(end+1)=TargetPresentationTime;
                    else
                        Target_Index{TargetType,1}(end+1)=TargetPresentationTime;
                    end
                end
            end            
            
            if strcmp(Stim_Time{1,1}(i,:), 'Attend2Audio_Fixation')
                
                InBlock = 0;
                BlockEnd = [];
                BlockType = [];
                CurrentCondition = 1;
                StimPresented = 0;
                BlockOrder(end+1,1:3) = [CurrentCondition 0 0];
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'Attend2Visual_Fixation')
                
                InBlock = 0;
                BlockEnd = [];
                BlockType = [];
                CurrentCondition = 2;
                StimPresented = 0;
                BlockOrder(end+1,1:3) = [CurrentCondition 0 0];
                
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'AudioOnly_Trial')
                if InBlock == 0
                    InBlock = 1;
                    BlockType = 1;
                    Block_Start{BlockType,CurrentCondition}(end+1) = Stim_Time{1,2}(i);
                    StimCount = StimCount + 1;
                    BlockOrder(end,2) = BlockType;
                else
                    StimCount = StimCount + 1;
                end
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'VisualOnly_Trial')
                if InBlock == 0
                    InBlock = 1;
                    BlockType = 2;
                    Block_Start{BlockType,CurrentCondition}(end+1) = Stim_Time{1,2}(i);
                    StimCount = StimCount + 1;
                    BlockOrder(end,2) = BlockType;
                else
                    StimCount = StimCount + 1;
                end
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'AudioVisual_Trial')
                if InBlock == 0
                    InBlock = 1;
                    BlockType = 3;
                    Block_Start{BlockType,CurrentCondition}(end+1) = Stim_Time{1,2}(i);
                    StimCount = StimCount + 1;
                    BlockOrder(end,2) = BlockType;
                else
                    StimCount = StimCount + 1;
                end
                
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'Auditory_Target')
                
                if StimPresented
                    if Target==0
                        Target_Index{TargetType,CurrentCondition}(end+1)=TargetPresentationTime;
                    else
                        if CurrentCondition==1
                            Target_Index{TargetType,2}(end+1)=TargetPresentationTime;
                        else
                            Target_Index{TargetType,1}(end+1)=TargetPresentationTime;
                        end
                    end
                end

                StimPresented = 1;
                TargetType = 1;
                if CurrentCondition == 1
                    Target = 1;
                else
                    Target = 0;
                end
                TargetPresentationTime = Stim_Time{1,2}(i);
                
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'Visual_Target')
                
                if StimPresented
                    if Target==0
                        Target_Index{TargetType,CurrentCondition}(end+1)=TargetPresentationTime;
                    else
                        if CurrentCondition==1
                            Target_Index{TargetType,2}(end+1)=TargetPresentationTime;
                        else
                            Target_Index{TargetType,1}(end+1)=TargetPresentationTime;
                        end
                    end
                end
                
                StimPresented = 1;
                TargetType = 2;
                if CurrentCondition == 2
                    Target = 1;
                else
                    Target = 0;
                end
                TargetPresentationTime = Stim_Time{1,2}(i);
                
            elseif strcmp(Stim_Time{1,1}(i,:), '1') && ~isempty(CurrentCondition)
                
                if BlockOrder(end,3)==0
                    BlockOrder(end,3)=1;
                end

                if StimPresented
                    Response_Index{TargetType,CurrentCondition}(end+1)=Stim_Time{1,2}(i);
                    if Target
                        Target_Index{TargetType,CurrentCondition}(end+1)=TargetPresentationTime;
                    else
                        if CurrentCondition==1
                            Target_Index{TargetType,2}(end+1)=TargetPresentationTime;
                        else
                            Target_Index{TargetType,1}(end+1)=TargetPresentationTime;
                        end
                    end
                        
                else
                    Extra_Response{1}(end+1)=Stim_Time{1,2}(i);
                end
                
                StimPresented = 0;
                
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'Stimulus_Offset')
                if StimCount == 50
                    BlockEnd = Stim_Time{1,2}(i);
                    Block_Duration{BlockType,CurrentCondition}(end+1) = ...
                        BlockEnd-Block_Start{BlockType,CurrentCondition}(end);
                    StimCount = 0;
                end
                
            elseif strcmp(Stim_Time{1,1}(i,:), 'Final_Fixation')
                if StimPresented
                    if Target==0
                        Target_Index{TargetType,CurrentCondition}(end+1)=TargetPresentationTime;
                    else
                        if CurrentCondition==1
                            Target_Index{TargetType,2}(end+1)=TargetPresentationTime;
                        else
                            Target_Index{TargetType,1}(end+1)=TargetPresentationTime;
                        end
                    end
                end                
                
                
            end
            
        end
        
        A = ~BlockOrder(:,3);
        B = [A [A(2:end); 0]];
        B = all(B,2);
        B(find(B)+1) = 1;
        BlockOrder(:,3)=B;

        for iStimCond=1:3
            for iAtt=1:2
                A = BlockOrder(all([BlockOrder(:,1)==iAtt,BlockOrder(:,2)==iStimCond],2),3);
                RespInBlock(iStimCond,iAtt,Block(iStimCond,iAtt):Block(iStimCond,iAtt)+2) = A;
                Block(iStimCond,iAtt)=Block(iStimCond,iAtt)+3;
                clear A
            end
        end    
        
        
        if sum(strcmp(Stim_Time{1,1}, 'Visual_Target')) + sum(strcmp(Stim_Time{1,1}, 'Auditory_Target')) ~= ...
                sum(sum(cellfun('length', Target_Index)))
            error('We are missing some targets')
        end
        
        
        %%
        %         Block_Start
        %         Block_Duration
        %
        %         Target_Index
        %         Response_Index
        
        %             'A Stim - Auditory Attention', ...
        %             'V Stim - Auditory Attention', ...
        %             'AV Stim - Auditory Attention', ...
        %             'A Stim - Visual Attention', ...
        %             'V Stim - Visual Attention', ...
        %             'AV Stim - Visual Attention', ...
        
        %             'A Targets - Auditory Attention', ...
        %             'A Targets - Visual Attention', ...
        %             'V Targets - Auditory Attention', ...
        %             'V Targets - Visual Attention', ...
        
        %             'Extra responses',...
        
        %             'Responses - Auditory Attention', ...
        %             'Responses - Visual Attention'
        
        figure('Name', ['Subject ' SubjID], 'Position', [100, 100, 1500, 1000])
%         subplot(NbRuns,1,RunInd)
        
        ylabel(['Run ' num2str(RunInd)])
        xlabel('Time (s)')
        
        hold on
        
        for i=1:numel(Block_Start)
            for j=1:numel(Block_Start{i})
                rectangle('Position', [Block_Start{i}(j) 0 Block_Duration{i}(j) .6], ...
                    'facecolor', COLOR(i,:), 'edgecolor', COLOR(i,:))
            end
        end
        
        ColorTarget = 'rgrg';
        Height = [0.7 0.8 0.7 0.8];
        for i=1:numel(Target_Index)
            stem(Target_Index{i},ones(size(Target_Index{i}))*Height(i), 'Color', ColorTarget(i),...
             'MarkerFaceColor', ColorTarget(i), 'MarkerSize', 4)
        end
        
        ColorResponse = 'bwwb';
        for i=1:numel(Response_Index)
            stem(Response_Index{i},ones(size(Response_Index{i})), 'Color', 'b',...
             'MarkerFaceColor', ColorResponse(i), 'MarkerSize', 4)
        end
        stem(Extra_Response{1},ones(size(Extra_Response{1}))*1.1, 'k', 'MarkerFaceColor', 'k', ...
                'MarkerSize', 4)

        A = axis;
        axis([A(1) A(2) A(3) 1.15])
        
        set(gca,'tickdir', 'out', 'xtick', 0:30:A(2), 'xticklabel',0:30:A(2),...
        'ytick', [], 'yticklabel', [], 'ticklength', [0.001 0.001], 'fontsize', 6)
    
        mtit(['Subject ' SubjID])
        
        Name = fullfile('/data','AV_Integration_2','Figures','Behavioral', ...
            ['ResultsTimingSubject' SubjID 'Run' num2str(RunInd) '.tif']);
        print(gcf, Name,'-dtiff')
        
    end

    close all
    
    ExcludeBlocks(SubjInd).ExcludeBlocks = RespInBlock;
   
%     mtit(['Subject ' SubjID])
    
%     Name = fullfile('/data','AV_Integration_2','Figures','Behavioral', ...
%         ['ResultsTimingSubject' SubjID '.tif']);
%     print(gcf, Name,'-dtiff')
    
    cd (StartFolder)
    
end

 
save(fullfile('/data','AV_Integration_2','Figures','Behavioral','ExcludeBlocks.mat'), 'ExcludeBlocks')

