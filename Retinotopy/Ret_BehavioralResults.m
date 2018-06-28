clc
clear all
close all

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

% Color for Conditions
COLOR =   [...
    255 150 150; ...
    150 255 150; ...
    150 220 255; ...
    255 75 75; ...
    75 255 75; ...
    75 75 255];
COLOR=COLOR/255;

% Color for Subjects
COLOR_Subject= [
    0,0,0;
    31,120,180;
    178,223,138;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

StartDirectory = pwd;

RT_HIT = cell(1,2,size(SubjectList,1));
HIT_TOTAL = nan(1,2,size(SubjectList,1));
MISS_TOTAL = nan(1,2,size(SubjectList,1));
EXTRA_ANSWERS_TOTAL = nan(1,2,size(SubjectList,1));

IndStart = 5;% first row of data points in txt file

TargetTimeOut = 2; % s
TargetTimeOut = TargetTimeOut * 10000;

ResponseTimeIn = 0.15;
ResponseTimeIn = ResponseTimeIn * 10000;

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\nRunning subject %s.\n\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    BehavioralFolder = fullfile(SubjectFolder, 'Retinotopy', 'Behavioral');
    
    cd(BehavioralFolder)
    
    LogFileList = dir('Logfile*.txt');
    
    %%
    for FileInd = 1:size(LogFileList,1)
        
        %%
        disp(LogFileList(FileInd).name)
        
        fid = fopen(fullfile (pwd, LogFileList(FileInd).name));
        FileContent = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s', 'headerlines', IndStart, 'returnOnError',0);
        fclose(fid);
        
        EOF = find(strcmp('Final_Fixation', FileContent{1,3}));
        if isempty(EOF)
            EOF = find(strcmp('Quit', FileContent{1,2})) - 1;
        end
        
        Stim_Time = {FileContent{1,3}(1:EOF)  FileContent{1,4}(1:EOF)};
        
        TEMP = [ ...
            find(strcmp('NewCycle', Stim_Time{1,1})) ;...
            find(strcmp('START', Stim_Time{1,1})) ;...
            find(strcmp('Instruction', Stim_Time{1,1})) ;...
            find(strcmp('30', Stim_Time{1,1})) ;...
            find(strcmp('Final_Fixation', Stim_Time{1,1}))];
        
        IndexTargets = find(strcmp('Visual_Target', Stim_Time{1,1}));
        
        Stim_Time{1,1}(TEMP,:) = [];
        Stim_Time{1,2}(TEMP,:) = [];
        
        %%
        HIT = 0;
        MISS = 0;
        EXTRA_ANSWERS = 0;
        
        StimPresented = 0;
        TargetPresentationTime = 0;
        
        % Loop to analyze the whole run;
        for Ind = 1 : length(Stim_Time{1,1})
            
            if StimPresented == 1 && str2double(char(Stim_Time{1,2}(Ind))) > TargetPresentationTime + TargetTimeOut
                StimPresented = 0;
                MISS = MISS + 1;
            end
            
            if strcmp('Visual_Target', Stim_Time{1,1}(Ind))
                if StimPresented == 1;
                    MISS = MISS + 1;
                end
                StimPresented = 1;
                TargetPresentationTime = str2double(char(Stim_Time{1,2}(Ind)));
                
            elseif strcmp('1', Stim_Time{1,1}(Ind))
                if str2double(char(Stim_Time{1,2}(Ind))) - TargetPresentationTime <  ResponseTimeIn
                    EXTRA_ANSWERS = EXTRA_ANSWERS + 1;
                elseif StimPresented == 1
                    HIT = HIT + 1;
                    RT_HIT{1, FileInd, SubjInd}(end+1) = (str2double(char(Stim_Time{1,2}(Ind)))-TargetPresentationTime)/10000;
                    StimPresented = 0;
                else
                    EXTRA_ANSWERS = EXTRA_ANSWERS + 1;
                end
            else
            end
            
        end
        
        
        if StimPresented == 1
            MISS = MISS + 1;
        end;
        
        if length(IndexTargets)~=HIT+MISS
            warning('Houston ! We are missing some targets !'); %#ok<WNTAG>
        end
        
        
        HIT_TOTAL(1,FileInd,SubjInd) = HIT;
        MISS_TOTAL(1,FileInd,SubjInd) = MISS;
        EXTRA_ANSWERS_TOTAL(1,FileInd,SubjInd) = EXTRA_ANSWERS;
        
    end
    
    fprintf(['\nSubject ' SubjectList(SubjInd,:) ' average\n'])
    
    fprintf(' Hits: %i\n', sum(HIT_TOTAL(1,:,SubjInd)))
    fprintf(' Misses: %i\n', sum(MISS_TOTAL(1,:,SubjInd)))
    fprintf(' Extra: %i\n\n', sum(EXTRA_ANSWERS_TOTAL(1,:,SubjInd)))
    
    
    cd(StartDirectory)
    
end

%% Plot perf
close all

LEGEND=char({'Hits';'Miss';'Extra'});

FontSize = 8;

figure('name', 'Retinotopy results', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

hold on
grid on

bar([squeeze(nansum(HIT_TOTAL,2)) ...
    squeeze(nansum(MISS_TOTAL,2)) ...
    squeeze(nansum(EXTRA_ANSWERS_TOTAL,2))])

set(gca, 'xtick',1:size(SubjectList,1), 'xticklabel', SubjectList, ...
    'fontsize',FontSize);

t=title('Retinotopy: behavioral performances.');
set(t,'fontsize',FontSize);

legend(LEGEND)

print(gcf, fullfile(StartDirectory, 'RetinotopyBehavioral.tif'), '-dtiff')



%% Reaction time

MAX = max([RT_HIT{:}])+.1;

figure('name', 'Reaction time', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])

hold on
grid on

for iSubj=1:size(SubjectList,1)
    
    tmp = [RT_HIT{1,1,iSubj} RT_HIT{1,2,iSubj}];
    
    Values2Plot{1,iSubj} = tmp;      
    
end

distributionPlot(Values2Plot, 'showMM', 1, 'histOpt', 1, ...
                                'globalNorm', 0, 'color', 'b')

set(gca, 'xtick',1:size(SubjectList,1), 'xticklabel', SubjectList, ...
    'fontsize',FontSize);

t=title('Retinotopy: RT.');
set(t,'fontsize',FontSize);

print(gcf, fullfile(StartDirectory,  'RetinotopyReactionTime.tif'), '-dtiff')

