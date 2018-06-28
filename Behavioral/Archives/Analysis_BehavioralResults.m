clc
clear
close all

SubjectList = [...
    '02';...
    '03';...
    '04';...
%     '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
%     '14';...
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
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

StartDirectory = fullfile(pwd, '..','..');
FiguresDirectory = fullfile(StartDirectory, 'Figures', 'Behavioral');

RT_HIT = cell(3,2,size(SubjectList,1));
RT_FA = cell(3,2,size(SubjectList,1));

IndStart = 5;% first row of data points in txt file

TargetTimeOut = 2.5; % s
TargetTimeOut = TargetTimeOut * 10000;

ResponseTimeIn = 0.15;
ResponseTimeIn = ResponseTimeIn * 10000;

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    fprintf('\nRunning subject %s.\n\n', SubjID)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    BehavioralFolder = fullfile(SubjectFolder, 'Runs');
    
    cd(BehavioralFolder)
    
    %%
    LogFileList = dir('Logfile*.txt');
    
    HIT_TOTAL = zeros(4,3,size(LogFileList,1));
    MISS_TOTAL = zeros(4,3,size(LogFileList,1));
    FALSE_ALARM_TOTAL = zeros(4,3,size(LogFileList,1));
    CORRECT_REJECTION_TOTAL = zeros(4,3,size(LogFileList,1));
    
    HIT_Block_TOTAL = [];
    MISS_Block_TOTAL = [];
    FALSE_ALARM_Block_TOTAL = [];
    CORRECT_REJECTION_Block_TOTAL = [];
        
    GroupResults(SubjInd).d = zeros(4,3,size(LogFileList,1));
    GroupResults(SubjInd).Accuracy = zeros(4,3,size(LogFileList,1));
    
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
            find(strcmp('SOA_Fix', Stim_Time{1,1})) ;...
            find(strcmp('SOA_Var', Stim_Time{1,1})) ;...
            find(strcmp('30', Stim_Time{1,1})) ;...
            find(strcmp('Final_Fixation', Stim_Time{1,1})) ;...
            find(strcmp('Fixation_Onset_Fix', Stim_Time{1,1})) ;...
            find(strcmp('Stimulus_Offset', Stim_Time{1,1})); ...
            find(strcmp('Long_Fixation', Stim_Time{1,1}))];
        
        IndexTargets = [...
            find(strcmp('Visual_Target', Stim_Time{1,1}));...
            find(strcmp('Auditory_Target', Stim_Time{1,1}))];
        
        
        Stim_Time{1,1}(TEMP,:) = [];
        Stim_Time{1,2}(TEMP,:) = [];
        
        %%
        HIT = zeros(4,3);
        MISS = zeros(4,3);
        FALSE_ALARM = zeros(4,3);
        CORRECT_REJECTION = zeros(4,3);
        
        HIT_Block = zeros(4,3,3);
        MISS_Block = zeros(4,3,3);
        FALSE_ALARM_Block = zeros(4,3,3);
        CORRECT_REJECTION_Block = zeros(4,3,3);
        
        EXTRA_ANSWERS = 0;
        
        StimPresented = 0;
        TargetORDistractor = 0;
        TargetPresentationTime = 0;
        
        Block =  zeros(3,2);
        
        % Loop to analyze the whole run;
        for Ind = 1 : length(Stim_Time{1,1})
            
            if StimPresented == 1 && str2double(char(Stim_Time{1,2}(Ind))) > TargetPresentationTime + TargetTimeOut
                StimPresented = 0;
                if TargetORDistractor == 1
                    MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
                    MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;
                else
                    CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
                    CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;
                end
            end
            
            
            if  strcmp('Attend2Audio_Fixation', Stim_Time{1,1}(Ind))
                CurrentCondition = 1; CurrentTrialType=0;
            elseif  strcmp( 'Attend2Visual_Fixation', Stim_Time{1,1}(Ind))
                CurrentCondition = 2; CurrentTrialType=0;
                
            elseif  strcmp('AudioOnly_Trial', Stim_Time{1,1}(Ind))
                if CurrentTrialType==0
                    Block(1, CurrentCondition) =  Block(1, CurrentCondition) + 1;
                end
                CurrentTrialType = 1;
            elseif strcmp('VisualOnly_Trial', Stim_Time{1,1}(Ind))
                if CurrentTrialType==0
                    Block(2, CurrentCondition) =  Block(2, CurrentCondition) + 1;
                end                
                CurrentTrialType = 2;
            elseif strcmp('AudioVisual_Trial', Stim_Time{1,1}(Ind))
                if CurrentTrialType==0
                    Block(3, CurrentCondition) =  Block(3, CurrentCondition) + 1;
                end                
                CurrentTrialType = 3;
            end

            
            if strcmp('Auditory_Target', Stim_Time{1,1}(Ind))
                if StimPresented == 1;
                    if TargetORDistractor == 1
                        MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
                    MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                        
                    else
                        CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
                    CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                         
                    end
                end
                StimPresented = 1; TrialTypeOfInterest = CurrentTrialType;
                TargetPresentationTime = str2double(char(Stim_Time{1,2}(Ind)));
                
                if CurrentCondition==1
                    TargetORDistractor = 1;
                else
                    TargetORDistractor = 0;
                end
                
            elseif strcmp('Visual_Target', Stim_Time{1,1}(Ind))
                if StimPresented == 1;
                    if TargetORDistractor == 1
                        MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
                    MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        MISS_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                                  
                    else
                        CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
                    CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        CORRECT_REJECTION_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                         
                    end
                end
                StimPresented = 1; TrialTypeOfInterest = CurrentTrialType;
                TargetPresentationTime = str2double(char(Stim_Time{1,2}(Ind)));
                
                if CurrentCondition==2
                    TargetORDistractor = 1;
                else
                    TargetORDistractor = 0;
                end
                
            elseif strcmp('1', Stim_Time{1,1}(Ind))
                if str2double(char(Stim_Time{1,2}(Ind))) - TargetPresentationTime <  ResponseTimeIn
                    EXTRA_ANSWERS = EXTRA_ANSWERS + 1;
                elseif StimPresented == 1
                    if TargetORDistractor == 1
                        HIT(TrialTypeOfInterest, CurrentCondition) = HIT(TrialTypeOfInterest, CurrentCondition) + 1;
                    HIT_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        HIT_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                                  
                        RT_HIT{TrialTypeOfInterest, CurrentCondition, SubjInd}(end+1) = (str2double(char(Stim_Time{1,2}(Ind)))-TargetPresentationTime)/10000;
                    else
                        FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) = FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) + 1;
                    FALSE_ALARM_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) = ...
                        FALSE_ALARM_Block(TrialTypeOfInterest, CurrentCondition, Block(TrialTypeOfInterest, CurrentCondition)) + 1;                           
                        RT_FA{TrialTypeOfInterest, CurrentCondition, SubjInd}(end+1) = (str2double(char(Stim_Time{1,2}(Ind)))-TargetPresentationTime)/10000;
                    end
                    StimPresented = 0;
                else
                    EXTRA_ANSWERS = EXTRA_ANSWERS + 1;
                end
            else
            end
            
        end
        
        ExtravResp(1,FileInd) = EXTRA_ANSWERS;
        
        if any([HIT_Block([1:3 5:7])+MISS_Block([1:3 5:7])]==0)
%             HIT_Block+MISS_Block
%             Stim_Time{1,1}
        end
        
        if StimPresented == 1
            MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
        end;
        
        if length(IndexTargets)~=sum(sum(HIT+MISS+FALSE_ALARM+CORRECT_REJECTION))
            warning('Houston ! We are missing some targets !'); %#ok<WNTAG>
        end
        
        %%
        for i=1:size(HIT,1)
            HIT(i,end) = sum(HIT(i,1:end-1));
            MISS(i,end) = sum(MISS(i,1:end-1));
            FALSE_ALARM(i,end) = sum(FALSE_ALARM(i,1:end-1));
            CORRECT_REJECTION(i,end) = sum(CORRECT_REJECTION(i,1:end-1));
        end
        
        for j=1:size(HIT,2)
            HIT(end,j) = sum(HIT(1:end-1,j));
            MISS(end,j) = sum(MISS(1:end-1,j));
            FALSE_ALARM(end,j) = sum(FALSE_ALARM(1:end-1,j));
            CORRECT_REJECTION(end,j) = sum(CORRECT_REJECTION(1:end-1,j));
        end
        
        HIT_TOTAL(:,:,FileInd) = HIT;
        MISS_TOTAL(:,:,FileInd) = MISS;
        FALSE_ALARM_TOTAL(:,:,FileInd) = FALSE_ALARM;
        CORRECT_REJECTION_TOTAL(:,:,FileInd) = CORRECT_REJECTION;
        
        HIT_Block_TOTAL = cat(3, HIT_Block_TOTAL, HIT_Block);
        MISS_Block_TOTAL = cat(3, MISS_Block_TOTAL, MISS_Block);
        FALSE_ALARM_Block_TOTAL = cat(3, FALSE_ALARM_Block_TOTAL, FALSE_ALARM_Block);
        CORRECT_REJECTION_Block_TOTAL = cat(3, CORRECT_REJECTION_Block_TOTAL, CORRECT_REJECTION_Block);
        
        
        for i=1:size(HIT,1)
            for j=1:size(HIT,2)
                
                FalseAlarmRate = FALSE_ALARM(i,j)/(FALSE_ALARM(i,j)+CORRECT_REJECTION(i,j));
                if FalseAlarmRate==1
                    FalseAlarmRate = 1 - 1/(2*(CORRECT_REJECTION(i,j)+FALSE_ALARM(i,j)));
                end
                if FalseAlarmRate==0
                    FalseAlarmRate = 1/(2*(CORRECT_REJECTION(i,j)+FALSE_ALARM(i,j)));
                end
                
                HitRate = HIT(i,j)/(HIT(i,j)+MISS(i,j));
                if HitRate==1
                    HitRate = 1 - 1/(2*((HIT(i,j)+MISS(i,j))));
                end
                if HitRate==0
                    HitRate = 1/(2*((HIT(i,j)+MISS(i,j))));
                end
                
                D_prime(i,j) = norminv(HitRate)-norminv(FalseAlarmRate);
                Accuracy(i,j) = (HIT(i,j) + CORRECT_REJECTION(i,j)) / (HIT(i,j)+MISS(i,j)+FALSE_ALARM(i,j)+CORRECT_REJECTION(i,j));
                
            end
        end
        
        GroupResults(SubjInd).d(:,:,FileInd) = D_prime;
        GroupResults(SubjInd).Accuracy(:,:,FileInd) = Accuracy;
        
        clear D_prime Accuracy
        
    end
    
    ExtravResp
    
    % Dirty trick in case "there was no target or distractor for a given
    % condition for a given block"
    HIT_Block_TOTAL((MISS_Block_TOTAL+HIT_Block_TOTAL)==0)=NaN;
    CORRECT_REJECTION_Block_TOTAL((FALSE_ALARM_Block_TOTAL+CORRECT_REJECTION_Block_TOTAL)==0)=NaN;
    
    fprintf('Subject %i average', SubjInd)
    
    GroupResults(SubjInd).HitRate(:,:,:) = HIT_TOTAL./(HIT_TOTAL+MISS_TOTAL); %#ok<*SAGROW>
    GroupResults(SubjInd).CorrectRejectionRate(:,:,:) = CORRECT_REJECTION_TOTAL./(CORRECT_REJECTION_TOTAL+FALSE_ALARM_TOTAL);
    GroupResults(SubjInd).HitRateBlock(:,:,:) = HIT_Block_TOTAL./(HIT_Block_TOTAL+MISS_Block_TOTAL); %#ok<*SAGROW>
    GroupResults(SubjInd).CorrectRejectionRateBlock(:,:,:) = CORRECT_REJECTION_Block_TOTAL./(CORRECT_REJECTION_Block_TOTAL+FALSE_ALARM_Block_TOTAL);
    
    HIT_TOTAL = sum(HIT_TOTAL,3);
    MISS_TOTAL = sum(MISS_TOTAL,3);
    FALSE_ALARM_TOTAL = sum(FALSE_ALARM_TOTAL,3);
    CORRECT_REJECTION_TOTAL = sum(CORRECT_REJECTION_TOTAL,3);
    
    for i=1:size(HIT_TOTAL,1)
        for j=1:size(HIT_TOTAL,2)
            
            FalseAlarmRate = FALSE_ALARM_TOTAL(i,j)/(FALSE_ALARM_TOTAL(i,j)+CORRECT_REJECTION_TOTAL(i,j));
            if FalseAlarmRate==1
                FalseAlarmRate = 1 - 1/(2*(CORRECT_REJECTION_TOTAL(i,j)+FALSE_ALARM_TOTAL(i,j)));
            end
            if FalseAlarmRate==0
                FalseAlarmRate = 1/(2*(CORRECT_REJECTION_TOTAL(i,j)+FALSE_ALARM_TOTAL(i,j)));
            end
            
            HitRate = HIT_TOTAL(i,j)/(HIT_TOTAL(i,j)+MISS_TOTAL(i,j));
            if HitRate==1
                HitRate = 1 - 1/(2*((HIT_TOTAL(i,j)+MISS_TOTAL(i,j))));
            end
            if HitRate==0
                HitRate = 1/(2*((HIT_TOTAL(i,j)+MISS_TOTAL(i,j))));
            end
            
            
            D_prime_TOTAL(i,j) = norminv(HitRate)-norminv(FalseAlarmRate);
            Accuracy_TOTAL(i,j) = (HIT_TOTAL(i,j) + CORRECT_REJECTION_TOTAL(i,j)) / (HIT_TOTAL(i,j)+MISS_TOTAL(i,j)+FALSE_ALARM_TOTAL(i,j)+CORRECT_REJECTION_TOTAL(i,j));
            
        end
    end
    
    D_prime_TOTAL
    Accuracy_TOTAL
    
    GroupResults(SubjInd).D_prime_TOTAL = D_prime_TOTAL;
    GroupResults(SubjInd).Accuracy_TOTAL = Accuracy_TOTAL;
    
    cd(StartDirectory)
    
end


FontSize = 6;

%% Accuracy
% 
% close all
% 
% TEMP = [];
% 
% 
% figure('name', 'Accuracy', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])
% 
% for iSubj=1:size(GroupResults,2)
%     TEMP(:,:,iSubj) = mean(GroupResults(iSubj).Accuracy,3);
%     
%     iSubplot = 1;
%     for iRow=1:3
%         for iCol=1:4
%             
%             subplot(3,4,iSubplot)
%             hold on
%             grid on
%             plot(repmat(iSubj, length(GroupResults(iSubj).Accuracy(iCol,iRow,:)), 1),...
%                 squeeze(GroupResults(iSubj).Accuracy(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
%                 'MarkerSize', 4)
%             iSubplot = iSubplot + 1;
%         end
%     end
%     
% end
% 
% TEMP2 = mean(TEMP,3);
% TEMP3 = std(TEMP,[],3);
% 
% iSubplot = 1;
% for iRow=1:3
%     for iCol=1:4
%         
%         subplot(3,4,iSubplot)
%         axis([0 size(GroupResults,2)+.5 0 1])
%         errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
%         set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:1, 'fontsize',FontSize-2);
%         iSubplot = iSubplot + 1;
%     end
% end
% 
% subplot(3,4,1)
% t=title('Auditory stimulation');
% set(t,'fontsize',FontSize);
% t=ylabel('Auditory attention');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,5)
% t=ylabel('Visual attention');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,9)
% t=ylabel('All');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,2)
% t=title('Visual stimulation');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,3)
% t=title('AudioVisual stimulation');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,4)
% t=title('All');
% set(t,'fontsize',FontSize);
% 
% print(gcf, fullfile(FiguresDirectory, 'Accuracy.tif'), '-dtiff')
% 
% clear TEMP TEMP2 TEMP3

%% D prime
% 
% TEMP = [];
% 
% figure('name', 'D prime', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])
% 
% for iSubj=1:size(GroupResults,2)
%     TEMP(:,:,iSubj) = mean(GroupResults(iSubj).d,3);
%     MIN(iSubj) = min(GroupResults(iSubj).d(:));
%     MAX(iSubj) = max(GroupResults(iSubj).d(:));
%     
%     iSubplot = 1;
%     for iRow=1:3
%         for iCol=1:4
%             
%             subplot(3,4,iSubplot)
%             hold on
%             grid on
%             plot(repmat(iSubj, length(GroupResults(iSubj).d(iCol,iRow,:)), 1), ...
%                 squeeze(GroupResults(iSubj).d(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
%                 'MarkerSize', 4)
%             iSubplot = iSubplot + 1;
%         end
%     end
%     
% end
% 
% TEMP2 = mean(TEMP,3);
% TEMP3 = std(TEMP,[],3);
% 
% MIN = min(MIN)-.1;
% MAX = max(MAX)+.1;
% 
% iSubplot = 1;
% for iRow=1:3
%     for iCol=1:4
%         
%         subplot(3,4,iSubplot)
%         plot([0 size(GroupResults,2)+.5],[0 0], '--k')
%         axis([0 size(GroupResults,2)+.5 MIN MAX])
%         errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
%         set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', -1:1:3,'fontsize',FontSize-2);
%         iSubplot = iSubplot + 1;
%     end
% end
% 
% subplot(3,4,1)
% t=title('Auditory stimulation');
% set(t,'fontsize',FontSize);
% t=ylabel('Auditory attention');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,5)
% t=ylabel('Visual attention');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,9)
% t=ylabel('All');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,2)
% t=title('Visual stimulation');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,3)
% t=title('AudioVisual stimulation');
% set(t,'fontsize',FontSize);
% 
% subplot(3,4,4)
% t=title('All');
% set(t,'fontsize',FontSize);
% 
% print(gcf, fullfile(FiguresDirectory,  'D_Prime.tif'), '-dtiff')



%% Reaction time

MAX = max([max([RT_HIT{:}]) max([RT_FA{:}])])+.1;

figure('name', 'Reaction time', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])

iSubplot = 1;
for iRow=1:2
    for iCol=1:3
        
        subplot(2,3,iSubplot)
        hold on
        grid on
        
        for iSubj=1:size(GroupResults,2)
            
            plot(repmat(iSubj-.2, numel(RT_HIT{iCol,iRow,iSubj}), 1), ...
                RT_HIT{iCol,iRow,iSubj}, ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 2)
            
            if ~isempty(RT_FA{iCol,iRow,iSubj})
                plot(repmat(iSubj+.2, numel(RT_FA{iCol,iRow,iSubj}), 1), ...
                    RT_FA{iCol,iRow,iSubj}, ' o', ...
                    'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                    'MarkerSize', 2)
            end
            
        end
        
        axis([0 size(GroupResults,2)+.5 0 MAX])
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:ceil(MAX),'fontsize',FontSize-2);
        
        iSubplot = iSubplot + 1;
    end
    
end

subplot(2,3,1)
t=title('Auditory stimulation');
set(t,'fontsize',FontSize);
t=ylabel('Auditory attention');
set(t,'fontsize',FontSize);

subplot(2,3,4)
t=ylabel('Visual attention');
set(t,'fontsize',FontSize);

subplot(2,3,2)
t=title('Visual stimulation');
set(t,'fontsize',FontSize);

subplot(2,3,3)
t=title('AudioVisual stimulation');
set(t,'fontsize',FontSize);


print(gcf, fullfile(FiguresDirectory,  'ReactionTime.tif'), '-dtiff')


%% Hit rate

TEMP = [];

figure('name', 'Hit rate', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = mean(GroupResults(iSubj).HitRate,3);
    MIN(iSubj) = min(GroupResults(iSubj).HitRate(:));
    MAX(iSubj) = max(GroupResults(iSubj).HitRate(:));
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            hold on
            grid on
            plot(repmat(iSubj, length(GroupResults(iSubj).HitRate(iCol,iRow,:)), 1), ...
                squeeze(GroupResults(iSubj).HitRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = mean(TEMP,3);
TEMP3 = std(TEMP,[],3);

MIN = min(MIN)-.1;
MAX = max(MAX)+.1;

iSubplot = 1;
for iRow=1:3
    for iCol=1:4
        
        subplot(3,4,iSubplot)
        plot([0 size(GroupResults,2)+.5],[0 0], '--k')
        axis([0 size(GroupResults,2)+.5 MIN MAX])
        errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:1,'fontsize',FontSize-2);
        iSubplot = iSubplot + 1;
    end
end

subplot(3,4,1)
t=title('Auditory stimulation');
set(t,'fontsize',FontSize);
t=ylabel('Auditory attention');
set(t,'fontsize',FontSize);

subplot(3,4,5)
t=ylabel('Visual attention');
set(t,'fontsize',FontSize);

subplot(3,4,9)
t=ylabel('All');
set(t,'fontsize',FontSize);

subplot(3,4,2)
t=title('Visual stimulation');
set(t,'fontsize',FontSize);

subplot(3,4,3)
t=title('AudioVisual stimulation');
set(t,'fontsize',FontSize);

subplot(3,4,4)
t=title('All');
set(t,'fontsize',FontSize);

print(gcf, fullfile(FiguresDirectory,  'HitRate.tif'), '-dtiff')


%% CR rate

TEMP = [];

figure('name', 'Correct Rejection Rate', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = mean(GroupResults(iSubj).CorrectRejectionRate,3);
    MIN(iSubj) = min(GroupResults(iSubj).CorrectRejectionRate(:));
    MAX(iSubj) = max(GroupResults(iSubj).CorrectRejectionRate(:));
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            hold on
            grid on
            plot(repmat(iSubj, length(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), 1), ...
                squeeze(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = mean(TEMP,3);
TEMP3 = std(TEMP,[],3);

MIN = min(MIN)-.1;
MAX = max(MAX)+.1;

iSubplot = 1;
for iRow=1:3
    for iCol=1:4
        
        subplot(3,4,iSubplot)
        plot([0 size(GroupResults,2)+.5],[0 0], '--k')
        axis([0 size(GroupResults,2)+.5 MIN MAX])
        errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:1,'fontsize',FontSize-2);
        iSubplot = iSubplot + 1;
    end
end

subplot(3,4,1)
t=title('Auditory stimulation');
set(t,'fontsize',FontSize);
t=ylabel('Auditory attention');
set(t,'fontsize',FontSize);

subplot(3,4,5)
t=ylabel('Visual attention');
set(t,'fontsize',FontSize);

subplot(3,4,9)
t=ylabel('All');
set(t,'fontsize',FontSize);

subplot(3,4,2)
t=title('Visual stimulation');
set(t,'fontsize',FontSize);

subplot(3,4,3)
t=title('AudioVisual stimulation');
set(t,'fontsize',FontSize);

subplot(3,4,4)
t=title('All');
set(t,'fontsize',FontSize);

print(gcf, fullfile(FiguresDirectory,  'CorrectRejectionRate.tif'), '-dtiff')



%% CR and hit rate per session/subject
% 
% close all
% 
% FontSize=10;
% 
% for iSubj=1:size(GroupResults,2)
%     
%     figure('name', ['Hit and CR Rate for subject ' SubjectList(iSubj,:)], ...
%         'position',  [100, 100, 1500, 1000], 'Color', [1 1 1], ...
%         'visible', 'off')
%     
%     iSubplot = 1;
%     for iRow=1:3
%         for iCol=1:4
%             
%             subplot(3,4,iSubplot)
%             
%             hold on
%             grid on
%             
%             
%             plot(1:size(GroupResults(iSubj).HitRate,3), squeeze(GroupResults(iSubj).HitRate(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
%                 'MarkerSize', 4)
%             
%             plot(.1+[1:size(GroupResults(iSubj).HitRate,3)], squeeze(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), ...
%                 'MarkerSize', 4)
%             
%             errorbar(0.35, mean(GroupResults(iSubj).HitRate(iCol,iRow,:)),...
%                 std(GroupResults(iSubj).HitRate(iCol,iRow,:)), ...
%                 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
% 
%             errorbar(0.45, mean(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)),...
%                 std(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ...
%                 'ob', 'MarkerSize', 5)
%             
%             if iSubplot==12
%                 legend({'Hit rate';'CR rate'; 'mean hit rate'; 'mean CR rate'}', ...
%                     'location', 'SouthEast', 'fontsize', FontSize-6)
%             end
%             
%             plot([0 size(GroupResults,2)+.5],[0.7 0.7], '--k')
%             
%             axis([0 size(GroupResults(iSubj).HitRate,3)+.5 0 1])
%             
%             set(gca, 'xtick',1:size(GroupResults(iSubj).HitRate,3), ...
%                 'xticklabel', 1:size(GroupResults(iSubj).HitRate,3), ...
%                 'ytick', 0:.2:1,'fontsize',FontSize);
%             
%             t=xlabel('Session');
%             set(t,'fontsize',FontSize);
%             
%             iSubplot = iSubplot + 1;
%         end
%     end
%         
%     
%     subplot(3,4,1)
%     t=title('Auditory stimulation');
%     set(t,'fontsize',FontSize-2);
%     t=ylabel('Auditory attention');
%     set(t,'fontsize',FontSize);
%     
%     subplot(3,4,5)
%     t=ylabel('Visual attention');
%     set(t,'fontsize',FontSize);
%     
%     subplot(3,4,9)
%     t=ylabel('All');
%     set(t,'fontsize',FontSize);
%     
%     subplot(3,4,2)
%     t=title('Visual stimulation');
%     set(t,'fontsize',FontSize-2);
%     
%     subplot(3,4,3)
%     t=title('AudioVisual stimulation');
%     set(t,'fontsize',FontSize-2);
%     
%     subplot(3,4,4)
%     t=title('All');
%     set(t,'fontsize',FontSize-2);
%     
%     print(gcf, fullfile(FiguresDirectory,...
%         ['Hit&CR_RateSubject_' SubjectList(iSubj,:) '.tif']), '-dtiff')
%     
%     
% end




%% CR and hit rate per session/subject
% 
% close all
% 
% FontSize=10;  
% 
% for iSubj=1:size(GroupResults,2)
%     
%     figure('name', ['Hit and CR Rate for subject ' SubjectList(iSubj,:)], 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])
%     
%     iSubplot = 1;
%     for iRow=1:2
%         for iCol=1:3
%             
%             subplot(2,3,iSubplot)
%             
%             hold on
%             grid on
%             
%             
%             plot(1:size(GroupResults(iSubj).HitRateBlock,3), ...
%                 squeeze(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
%                 'MarkerSize', 4)
%             
%             plot(.2+[1:size(GroupResults(iSubj).HitRateBlock,3)], ...
%                 squeeze(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)), ' o', ...
%                 'Color', COLOR_Subject(iSubj,:), ...
%                 'MarkerSize', 4)
%             
%             errorbar(0.35, mean(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)),...
%                 std(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)), ...
%                 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
% 
%             errorbar(0.45, mean(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)),...
%                 std(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)), ...
%                 'ob', 'MarkerSize', 5)
%             
%             if iSubplot==6
%                 legend({'Hit rate';'CR rate'; 'mean hit rate'; 'mean CR rate'}', ...
%                     'location', 'SouthEast', 'fontsize', FontSize-6)
%             end
%             
%             plot([0 size(GroupResults(iSubj).HitRateBlock,3)+.5],[0.7 0.7], '--k')
%             plot([3.5 3.5],[0 1], '-k', 'linewidth', 2)
%             plot([6.5 6.5],[0 1], '-k', 'linewidth', 2)
%             plot([9.5 9.5],[0 1], '-k', 'linewidth', 2)
%             
%             axis([0 size(GroupResults(iSubj).HitRateBlock,3)+.5 0 1])
%             
%             set(gca, 'xtick',1:size(GroupResults(iSubj).HitRateBlock,3), ...
%                 'xticklabel', 1:size(GroupResults(iSubj).HitRateBlock,3), ...
%                 'ytick', 0:.2:1,'fontsize',FontSize);
%             
%             t=xlabel('Blocks');
%             set(t,'fontsize',FontSize);
%             
%             iSubplot = iSubplot + 1;
%         end
%     end
%         
%     
%     subplot(2,3,1)
%     t=title('Auditory stimulation');
%     set(t,'fontsize',FontSize-2);
%     t=ylabel('Auditory attention');
%     set(t,'fontsize',FontSize);
%     
%     subplot(2,3,4)
%     t=ylabel('Visual attention');
%     set(t,'fontsize',FontSize);
%     
%     subplot(2,3,2)
%     t=title('Visual stimulation');
%     set(t,'fontsize',FontSize-2);
%     
%     subplot(2,3,3)
%     t=title('AudioVisual stimulation');
%     set(t,'fontsize',FontSize-2);
%     
%     print(gcf, fullfile(FiguresDirectory,...
%         ['Hit&CR_Rate_PerBlock_Subject_' SubjectList(iSubj,:) '.tif']), '-dtiff')
%     
%     
% end
% 
