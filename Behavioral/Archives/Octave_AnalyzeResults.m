clc
clear all
close all

NbSubject = 16;

D_prime_GROUP = nan(4,3,length(NbSubject));
Accuracy_GROUP = nan(4,3,length(NbSubject));

RT_HIT_GROUP = nan(4,3,length(NbSubject));
RT_FALSE_ALARM_GROUP = nan(4,3,length(NbSubject));

StartDirectory = pwd;

for SubjInd = 1 : length(NbSubject)
    
    cd(strcat('Subject_', num2str(NbSubject(SubjInd))))
    
    cd('Runs')

    %%
    LogFileList = dir('Logfile*.txt');
    
    HIT_TOTAL = zeros(4,3,size(LogFileList,1));
    MISS_TOTAL = zeros(4,3,size(LogFileList,1));
    FALSE_ALARM_TOTAL = zeros(4,3,size(LogFileList,1));
    CORRECT_REJECTION_TOTAL = zeros(4,3,size(LogFileList,1));
    
    EXTRA_ANSWERS_TOTAL = 0;
    
    RT_HIT = cell(4,3);
    RT_FALSE_ALARM = cell(4,3);

    IndStart = 5;% first row of data points in txt file

    TargetTimeOut = 2.5; % in secs
    TargetTimeOut = TargetTimeOut * 10000;

    %%
    for FileInd = 1:size(LogFileList,1)

        %%
        disp(LogFileList(FileInd).name)

        fid = fopen(fullfile (pwd, LogFileList(FileInd).name));
        
        %octave
        FileContent = textscan(fid, '%s', 'delimiter','\t', 'headerlines', IndStart, 'returnOnError', 0);
        
        fclose(fid);
        
        
        EOF = find(strcmp('Final_Fixation', FileContent{1,1}));
        if isempty(EOF)
            EOF = find(strcmp('Quit', FileContent{1,1}));
        end

        Stim_Time = FileContent{1,1}(1:EOF(1)-3);


        
        % For Octave
        IndexStartAudioAttention = find(strcmp('Attend2Audio_Fixation', Stim_Time));
        
        IndexStartVisualAttention = find(strcmp('Attend2Visual_Fixation', Stim_Time));

        TEMP = find(strcmp('AudioOnly_Trial', Stim_Time));
        ToKeep = [IndexStartAudioAttention ; IndexStartVisualAttention ; TEMP(:)];
        clear TEMP

        TEMP = find(strcmp('VisualOnly_Trial', Stim_Time));
        ToKeep = [ToKeep ; TEMP(:)];
        clear TEMP

        TEMP = find(strcmp('AudioVisual_Trial', Stim_Time));
        ToKeep = [ToKeep ; TEMP(:)];
        clear TEMP

        TEMP = find(strcmp('Visual_Target', Stim_Time));
        ToKeep = [ToKeep ; TEMP(:)];
        clear TEMP

        TEMP = find(strcmp('Auditory_Target', Stim_Time));
        ToKeep = [ToKeep ; TEMP(:)];
        clear TEMP
        
        IndexResponses = find(strcmp('Response', Stim_Time));
        ToKeep = [ToKeep ; IndexResponses(:)+1];
        clear TEMP    

        Stim_Time1 = {Stim_Time{sort(ToKeep(:))}}';
    
        Stim_Time2 = {Stim_Time{sort(ToKeep(:))+1}}';
        
        clear ToKeep


	        IndexTargets = [...
            find(strcmp('Visual_Target', Stim_Time{1,1}));...
            find(strcmp('Auditory_Target', Stim_Time{1,1}))];
        
        Stim_Time{1,1} = Stim_Time1;
        Stim_Time{1,2} = Stim_Time2;

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
                        RT = (str2double(char(Stim_Time{1,2}(Ind))) - TargetPresentationTime)/10;
                    if TargetORDistractor == 1
                        HIT(TrialTypeOfInterest, CurrentCondition) = HIT(TrialTypeOfInterest, CurrentCondition) + 1;
                        RT_HIT{TrialTypeOfInterest,CurrentCondition}(end+1) = RT;
                        RT_HIT{TrialTypeOfInterest,3}(end+1) = RT;
                        RT_HIT{4,CurrentCondition}(end+1) = RT;
                        RT_HIT{4,3}(end+1) = RT;
                    else
                        FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) = FALSE_ALARM(TrialTypeOfInterest, CurrentCondition) + 1;
                        RT_FALSE_ALARM{TrialTypeOfInterest,CurrentCondition}(end+1) = RT;
                        RT_FALSE_ALARM{TrialTypeOfInterest,3}(end+1) = RT;
                        RT_FALSE_ALARM{4,CurrentCondition}(end+1) = RT;
                        RT_FALSE_ALARM{4,3}(end+1) = RT;
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
        
        EXTRA_ANSWERS_TOTAL = EXTRA_ANSWERS_TOTAL + EXTRA_ANSWERS;
        
        HIT;
        MISS;
        FALSE_ALARM;
        CORRECT_REJECTION;

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

        D_prime
        Accuracy
        
    end
    
    fprintf('Subject %i average\n', NbSubject(SubjInd))
    
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
    EXTRA_ANSWERS_TOTAL
    
    D_prime_GROUP(:,:,SubjInd) = D_prime_TOTAL;
    Accuracy_GROUP(:,:,SubjInd) = Accuracy_TOTAL;
    
    for i=1:size(RT_HIT,1)
        for j=1:size(RT_HIT,2)
                if isempty(RT_HIT{i,j})
                        TEMP1(i,j) = NaN;
                else
                        TEMP1(i,j) = median(RT_HIT{i,j});
                end
                
                if isempty(RT_FALSE_ALARM{i,j})
                        TEMP2(i,j) = NaN;
                else
                        TEMP2(i,j) = median(RT_FALSE_ALARM{i,j});
                end
        end
    end

    RT_HIT_GROUP(:,:,SubjInd) = TEMP1;
    RT_FALSE_ALARM_GROUP(:,:,SubjInd) = TEMP2;
    
    clear TEMP1 TEMP2
                                                   
    %figure(SubjInd)
    %subplot(221)
    %hist(RT_HIT{1,1})
    %subplot(222)
    %hist(RT_HIT{1,2})
    %subplot(223)
    %hist(RT_FALSE_ALARM{1,1})
    %subplot(224)
    %hist(RT_FALSE_ALARM{1,2})
    
    cd(StartDirectory)

end

return

%%
Legend = {'Audio ' ; 'Visual' ; 'Total '};

figure('name', 'd prime', 'color', 'w', 'position',  [1 1 500 800])
for i=1:4
        subplot(4,1,i)
        hold on
        box off
        grid on
        [H,P,CI(1,:)] = ttest(D_prime_GROUP(i,1,:));
        [H,P,CI(2,:)] = ttest(D_prime_GROUP(i,2,:));
        [H,P,CI(3,:)] = ttest(D_prime_GROUP(i,3,:));
        %errorbar(1:3, nanmean(D_prime_GROUP(i,:,:), 3), CI(:,1), CI(:,2), 'o ');
        plot(1:3, nanmean(D_prime_GROUP(i,:,:), 3), 'o ', 'MarkerFaceColor',[0 0 1]);
        plot([1 1], CI(1,:), '.-');
        plot([2 2], CI(2,:), '.-');
        plot([3 3], CI(3,:), '.-');
        scatter(ones(1, length(D_prime_GROUP(i,1,:)))+0.1, D_prime_GROUP(i,1,:), 'b');
        scatter(ones(1, length(D_prime_GROUP(i,2,:)))*2+0.1, D_prime_GROUP(i,2,:), 'b');
        scatter(ones(1, length(D_prime_GROUP(i,3,:)))*3+0.1, D_prime_GROUP(i,3,:), 'b');
        axis([0.4 3.5 0 ceil(max(max(max(D_prime_GROUP))))]);
        set(gca, 'xtick', 1:3 ,'xticklabel', [],  ...
                 'ytick', 0:ceil(max(max(max(D_prime_GROUP)))) ,'yticklabel', 0:ceil(max(max(max(D_prime_GROUP)))), ...
                 'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off') ;
end

subplot(4,1,4)
set(gca, 'xtick', 1:3 ,'xticklabel', Legend,  ...
         'ytick', 0:ceil(max(max(max(D_prime_GROUP)))) ,'yticklabel', 0:ceil(max(max(max(D_prime_GROUP)))), ...
         'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off') ;


ylabel('All trials', 'fontsize', 14);
xlabel('Attention Condition', 'fontsize', 14);

subplot(4,1,1)
ylabel('A trials', 'fontsize', 14);

subplot(4,1,2)
ylabel('V trials', 'fontsize', 14);

subplot(4,1,3)
ylabel('AV trials', 'fontsize', 14);      
     
print(gcf, 'd prime.tif', '-dtiff')     

%%
figure('name', 'accuracy', 'color', 'w', 'position',  [1 1 500 800])
for i=1:4
        subplot(4,1,i)
        hold on
        box off
        grid on
        [H,P,CI(1,:)] = ttest(Accuracy_GROUP(i,1,:));
        [H,P,CI(2,:)] = ttest(Accuracy_GROUP(i,2,:));
        [H,P,CI(3,:)] = ttest(Accuracy_GROUP(i,3,:));
        plot(1:3, nanmean(Accuracy_GROUP(i,:,:), 3), 'o ', 'MarkerFaceColor',[0 0 1]);
        plot([1 1], CI(1,:), '.-');
        plot([2 2], CI(2,:), '.-');
        plot([3 3], CI(3,:), '.-');
        %errorbar(1:3, nanmean(Accuracy_GROUP(i,:,:), 3), nanstd(Accuracy_GROUP(i,:,:), 3, 0), 'o ');
        scatter(ones(1, length(Accuracy_GROUP(i,1,:)))+0.1, Accuracy_GROUP(i,1,:), 'b');
        scatter(ones(1, length(Accuracy_GROUP(i,2,:)))*2+0.1, Accuracy_GROUP(i,2,:), 'b');
        scatter(ones(1, length(Accuracy_GROUP(i,3,:)))*3+0.1, Accuracy_GROUP(i,3,:), 'b');
        axis([0.4 3.5 0 1]);
        set(gca, 'xtick', 1:3 ,'xticklabel', [],  ...
                 'ytick', 0:0.2:1 ,'yticklabel', 0:0.2:1, ...
                 'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off') ;
end

subplot(4,1,4)
set(gca, 'xtick', 1:3 ,'xticklabel', Legend,  ...
         'ytick', 0:0.2:1 ,'yticklabel', 0:0.2:1, ...
         'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off');
     
ylabel('All trials', 'fontsize', 14);
xlabel('Attention Condition', 'fontsize', 14);

subplot(4,1,1)
ylabel('A trials', 'fontsize', 14);

subplot(4,1,2)
ylabel('V trials', 'fontsize', 14);

subplot(4,1,3)
ylabel('AV trials', 'fontsize', 14);
  
print(gcf, 'accuracy.tif', '-dtiff')

%%
figure('name', 'Reaction times', 'color', 'w', 'position',  [1 1 500 800])
for i=1:4
        subplot(4,1,i)
        hold on
        box off
        grid on
        
        plot([1:3]-0.1, nanmean(RT_HIT_GROUP(i,:,:), 3), 'o g', 'MarkerFaceColor',[0 1 0]);
        plot([1:3]+0.1, nanmean(RT_FALSE_ALARM_GROUP(i,:,:), 3), 'o r', 'MarkerFaceColor',[1 0 0]);
        
        [H,P,CI(1,:)] = ttest(RT_HIT_GROUP(i,1,:));
        [H,P,CI(2,:)] = ttest(RT_HIT_GROUP(i,2,:));
        [H,P,CI(3,:)] = ttest(RT_HIT_GROUP(i,3,:));
        plot([1 1]-0.1, CI(1,:), '.-g');
        plot([2 2]-0.1, CI(2,:), '.-g');
        plot([3 3]-0.1, CI(3,:), '.-g');
        %errorbar([1:3]-0.1, nanmean(RT_HIT_GROUP(i,:,:), 3), nanstd(RT_HIT_GROUP(i,:,:), 3, 0), 'o g');
        scatter(ones(1, length(RT_HIT_GROUP(i,1,:)))-0.05, RT_HIT_GROUP(i,1,:), 'g');
        scatter(ones(1, length(RT_HIT_GROUP(i,2,:)))*2-0.05, RT_HIT_GROUP(i,2,:), 'g');
        scatter(ones(1, length(RT_HIT_GROUP(i,3,:)))*3-0.05, RT_HIT_GROUP(i,3,:), 'g');
        
        [H,P,CI(1,:)] = ttest(RT_FALSE_ALARM_GROUP(i,1,:));
        [H,P,CI(2,:)] = ttest(RT_FALSE_ALARM_GROUP(i,2,:));
        [H,P,CI(3,:)] = ttest(RT_FALSE_ALARM_GROUP(i,3,:));
        plot([1 1]+0.1, CI(1,:), '.-r');
        plot([2 2]+0.1, CI(2,:), '.-r');
        plot([3 3]+0.1, CI(3,:), '.-r');
        %errorbar([1:3]+0.1, nanmean(RT_FALSE_ALARM_GROUP(i,:,:), 3), nanstd(RT_FALSE_ALARM_GROUP(i,:,:), 3, 0), 'o r');
        scatter(ones(1, length(RT_FALSE_ALARM_GROUP(i,1,:)))+0.15, RT_FALSE_ALARM_GROUP(i,1,:), 'r');
        scatter(ones(1, length(RT_FALSE_ALARM_GROUP(i,2,:)))*2+0.15, RT_FALSE_ALARM_GROUP(i,2,:), 'r');
        scatter(ones(1, length(RT_FALSE_ALARM_GROUP(i,3,:)))*3+0.15, RT_FALSE_ALARM_GROUP(i,3,:), 'r');
        axis([0.5 3.5 0 max(max(max([RT_HIT_GROUP RT_FALSE_ALARM_GROUP])))]);
        set(gca, 'xtick', 1:3 ,'xticklabel', [],  ...
                 'ytick', 0:250:1500 ,'yticklabel', 0:250:1500, ...
                 'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off') ;
end

subplot(4,1,4)
set(gca, 'xtick', 1:3 ,'xticklabel', Legend,  ...
         'ytick', 0:250:2000 ,'yticklabel', 0:250:2000, ...
         'ticklength', [0.01 0.01], 'tickdir', 'out', 'fontsize', 12, 'XGrid','off')
     
ylabel('All trials', 'fontsize', 14);
xlabel('Attention Condition', 'fontsize', 14);

legend(char({'Hits', 'False alarms'}), 'Location','NorthEast')

subplot(4,1,1)
ylabel('A trials', 'fontsize', 14);

subplot(4,1,2)
ylabel('V trials', 'fontsize', 14);

subplot(4,1,3)
ylabel('AV trials', 'fontsize', 14);  

print(gcf, 'RT.tif', '-dtiff')