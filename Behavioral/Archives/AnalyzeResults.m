clc
clear
close all

Subjects = [2 3 4 6 7 8 9 11 12 13 14 15 16];

StartDirectory = pwd;

for SubjInd = 1:length(Subjects)
    cd('Behavioral')

    cd(strcat('Subject_', num2str(Subjects(SubjInd))))

    cd('Runs')

    %%
    LogFileList = dir('Logfile*.txt');

    HIT_TOTAL = zeros(4,3,size(LogFileList,1));
    MISS_TOTAL = zeros(4,3,size(LogFileList,1));
    FALSE_ALARM_TOTAL = zeros(4,3,size(LogFileList,1));
    CORRECT_REJECTION_TOTAL = zeros(4,3,size(LogFileList,1));

    IndStart = 5;% first row of data points in txt file

    TargetTimeOut = 1000; % ms
    TargetTimeOut = TargetTimeOut * 10000;

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
        EXTRA_ANSWERS = 0;

        StimPresented = 0;
        TargetORDistractor = 0;
        TargetPresentationTime = 0;

        % Loop to analyze the whole run;
        for Ind = 1 : length(Stim_Time{1,1})

            if  strcmp('Attend2Audio_Fixation', Stim_Time{1,1}(Ind))
                CurrentCondition = 1;
            elseif  strcmp( 'Attend2Visual_Fixation', Stim_Time{1,1}(Ind))
                CurrentCondition = 2;

            elseif  strcmp('AudioOnly_Trial', Stim_Time{1,1}(Ind))
                CurrentTrialType = 1;
            elseif strcmp('VisualOnly_Trial', Stim_Time{1,1}(Ind))
                CurrentTrialType = 2;
            elseif strcmp('AudioVisual_Trial', Stim_Time{1,1}(Ind))
                CurrentTrialType = 3;
            end

            if StimPresented == 1 && str2double(char(Stim_Time{1,2}(Ind))) > TargetPresentationTime + TargetTimeOut
                StimPresented = 0;
                if TargetORDistractor == 1
                    MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
                else
                    CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
                end
            end

            if strcmp('Auditory_Target', Stim_Time{1,1}(Ind))
                if StimPresented == 1;
                    if TargetORDistractor == 1
                        MISS(TrialTypeOfInterest, CurrentCondition) = MISS(TrialTypeOfInterest, CurrentCondition) + 1;
                    else
                        CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
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
                    else
                        CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) = CORRECT_REJECTION(TrialTypeOfInterest, CurrentCondition) + 1;
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
                if StimPresented == 1
                    if TargetORDistractor == 1
                        HIT(TrialTypeOfInterest, CurrentCondition) = HIT(TrialTypeOfInterest, CurrentCondition) + 1;
                    else
                        FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) = FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) + 1;
                    end
                else
                    EXTRA_ANSWERS = EXTRA_ANSWERS + 1;
                end
                StimPresented = 0;
            else
            end

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

        %         D_prime
        %         Accuracy

        GroupResults(SubjInd).d(:,:,FileInd) = D_prime;
        GroupResults(SubjInd).Accuracy(:,:,FileInd) = Accuracy;

    end

    fprintf('Subject %i average', SubjInd)
    
    GroupResults(SubjInd).HitRate(:,:,:) = HIT_TOTAL./(HIT_TOTAL+MISS_TOTAL);
    GroupResults(SubjInd).CorrectRejectionRate(:,:,:) = CORRECT_REJECTION_TOTAL./(CORRECT_REJECTION_TOTAL+FALSE_ALARM_TOTAL);

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

%% Accuracy

TEMP = [];

figure('name', 'Accuracy', 'position', [100 100 800 500])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = GroupResults(iSubj).Accuracy_TOTAL;

    for iCol=1:3
        for iRow=1:4
            subplot(4,3,iCol+(iRow-1)*3)
            hold on
            grid on
            plot(repmat(iSubj, length(GroupResults(iSubj).Accuracy(iRow,iCol,:)), 1), squeeze(GroupResults(iSubj).Accuracy(iRow,iCol,:)), ' .')
        end
    end

end

TEMP2 = mean(TEMP,3)

for iCol=1:3
    for iRow=1:4
        subplot(4,3,iCol+(iRow-1)*3)
        axis([0 8 0 1])
        plot(0.5, TEMP2(iRow,iCol), 'ok')
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', Subjects, 'ytick', 0:.2:1);
    end
end


%% D prime

TEMP = [];

figure('name', 'D prime', 'position', [100 100 800 500])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = GroupResults(iSubj).D_prime_TOTAL;

    for iCol=1:3
        for iRow=1:4
            subplot(4,3,iCol+(iRow-1)*3)
            hold on
            grid on
            plot(repmat(iSubj, length(GroupResults(iSubj).d(iRow,iCol,:)), 1), squeeze(GroupResults(iSubj).d(iRow,iCol,:)), ' .')
        end
    end

end

TEMP2 = mean(TEMP,3)

for iCol=1:3
    for iRow=1:4
        subplot(4,3,iCol+(iRow-1)*3)
        axis([0 8 -2 4])
        plot(0.5, TEMP2(iRow,iCol), 'ok')
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', Subjects, 'ytick', -2:1:3);
    end
end


%% Hit Rate

TEMP = [];

figure('name', 'hit rate', 'position', [100 100 800 500])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = GroupResults(iSubj).HitRate;

    for iCol=1:3
        for iRow=1:4
            subplot(4,3,iCol+(iRow-1)*3)
            hold on
            grid on
            plot(repmat(iSubj, length(GroupResults(iSubj).HitRate(iRow,iCol,:)), 1), squeeze(GroupResults(iSubj).HitRate(iRow,iCol,:)), ' .')
        end
    end

end

TEMP2 = mean(TEMP,3)

for iCol=1:3
    for iRow=1:4
        subplot(4,3,iCol+(iRow-1)*3)
        axis([0 8 0 1])
        plot(0.5, TEMP2(iRow,iCol), 'ok')
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', Subjects, 'ytick', 0:.2:1);
    end
end