clc
clear
close all

StartDirectory = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

CodeFodler = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFodler, 'SubFun')))
Get_dependencies('/home/remi/Dropbox')

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

% The following lines are manually dirty hack to compensate for
% Presentation not logging the display of all targets
Correction_CR = cell(size(SubjectList,1),4);
Correction_CR{1,1} = [2,2,2,1];
Correction_CR{1,3} = [2,1,2,1];
Correction_CR{3,3} = [2,1,3,1];
Correction_CR{5,2} = [3,2,1,1;...
    2,1,1,1];
Correction_CR{5,3} = [2,2,3,1];
Correction_CR{6,1} = [3,2,3,1];
Correction_CR{6,2} = [3,1,3,1];
Correction_CR{6,3} = [2,1,2,1];
Correction_CR{6,4} = [2,2,2,1];
Correction_CR{7,3} = [3,1,1,1];
Correction_CR{7,4} = [3,2,1,1];
Correction_CR{8,2} = [1,1,2,1;...
    3,1,2,1];
Correction_CR{8,3} = [1,1,2,1];
Correction_CR{9,1} = [1,1,3,1];
Correction_CR{9,4} = [2,2,3,1];
Correction_CR{10,1} = [1,1,1,1];
Correction_CR{10,2} = [3,2,1,1;...
    1,1,2,1];
Correction_CR{12,1} = [1,1,1,1;...
    1,2,2,1];
Correction_CR{12,2} = [1,2,2,1];
Correction_CR{13,2} = [1,1,2,1];
Correction_CR{13,3} = [2,1,3,1];
Correction_CR{13,4} = [1,2,2,1];

Correction_Miss = cell(size(SubjectList,1),4);
Correction_Miss{1,3} = [2,2,3,1];
Correction_Miss{2,2} = [3,2,3,1];
Correction_Miss{2,4} = [2,2,2,1];
Correction_Miss{5,3} = [3,2,2,1];
Correction_Miss{7,4} = [1,2,1,1];
Correction_Miss{8,1} = [1,2,1,1];
Correction_Miss{13,2} = [1,2,1,1];

Correction_Hit = cell(size(SubjectList,1),4);
Correction_Hit{1,2} = [1,1,3,2];
Correction_Hit{6,2} = [2,1,1,1];
Correction_Hit{6,3} = [3,1,2,1];
Correction_Hit{7,4} = [3,1,3,1];
Correction_Hit{9,2} = [2,1,1,1];
Correction_Hit{10,1} = [2,1,2,1];
Correction_Hit{12,2} = [1,1,1,1];
Correction_Hit{12,3} = [1,1,1,1;...
    3,1,2,1];
Correction_Hit{13,1} = [1,1,3,1];

Correction_FA = cell(size(SubjectList,1),4);
Correction_FA{13,3} = [3,2,2,1];




% Color for Conditions
COLOR =   [...
    255 150 150; ...
    150 255 150; ...
    150 220 255; ...
    255 75 75; ...
    75 255 75; ...
    75 75 255];
COLOR=COLOR/255;

COLOR_Subject = ColorSubject();

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
    BehavioralFolder = fullfile(SubjectFolder, 'Behavioral', 'Runs');
    
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
        
        if StimPresented == 1
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
        
        %Correction_CR{1,1} = [2,2,2]
        
        if ~isempty(Correction_CR{SubjInd,FileInd})
            for i=1:size(Correction_CR{SubjInd,FileInd},1)
                CORRECT_REJECTION_Block(Correction_CR{SubjInd,FileInd}(i,1),...
                    Correction_CR{SubjInd,FileInd}(i,2),...
                    Correction_CR{SubjInd,FileInd}(i,3)) = Correction_CR{SubjInd,FileInd}(i,4);
            end
        end
        
        if ~isempty(Correction_Miss{SubjInd,FileInd})
            MISS_Block(Correction_Miss{SubjInd,FileInd}(1),...
                Correction_Miss{SubjInd,FileInd}(2),...
                Correction_Miss{SubjInd,FileInd}(3)) = Correction_Miss{SubjInd,FileInd}(4);
        end
        
        if ~isempty(Correction_FA{SubjInd,FileInd})
            FALSE_ALARM_Block(Correction_FA{SubjInd,FileInd}(1),...
                Correction_FA{SubjInd,FileInd}(2),...
                Correction_FA{SubjInd,FileInd}(3)) = Correction_FA{SubjInd,FileInd}(4);
        end
        
        if ~isempty(Correction_Hit{SubjInd,FileInd})
            for i=1:size(Correction_Hit{SubjInd,FileInd},1)
                HIT_Block(Correction_Hit{SubjInd,FileInd}(i,1),...
                    Correction_Hit{SubjInd,FileInd}(i,2),...
                    Correction_Hit{SubjInd,FileInd}(i,3)) = Correction_Hit{SubjInd,FileInd}(i,4);
            end
        end
        
        A = [HIT_Block(1:3,1:2,:)+MISS_Block(1:3,1:2,:)]==0;
        B = [CORRECT_REJECTION_Block(1:3,1:2,:)+FALSE_ALARM_Block(1:3,1:2,:)]==0;
        
        
        
        if any(A(:)) || any(B(:))
            
            fprintf('\nSubject: %s ;  Run: %i \n\n', SubjID, FileInd)
            
            if any(A(:))
                tmp = zeros(size(HIT_Block));
                tmp(1:3,1:2,:) = A;
                
                HIT_Block(logical(tmp)) = NaN;
                MISS_Block(logical(tmp)) = NaN;
                A;
            end
            
            if any(B(:))
                tmp = zeros(size(HIT_Block));
                tmp(1:3,1:2,:) = B;
                
                FALSE_ALARM_Block(logical(tmp)) = NaN;
                CORRECT_REJECTION_Block(logical(tmp)) = NaN;
                B;
            end
            
            
            Stim_Time{1,1};
            clear A B
            
            
        end
        clear A B
        
        
        if length(IndexTargets)~=sum(sum(HIT+MISS+FALSE_ALARM+CORRECT_REJECTION))
            warning('Houston ! We are missing some targets !'); %#ok<WNTAG>
        end
        
        HIT = nansum(HIT_Block,3);
        MISS = nansum(MISS_Block,3);
        FALSE_ALARM = nansum(FALSE_ALARM_Block,3);
        CORRECT_REJECTION = nansum(CORRECT_REJECTION_Block,3);
        
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
    
    GroupResults(SubjInd).RT_FA = RT_FA(:,:,SubjInd);
    GroupResults(SubjInd).RT_HIT = RT_HIT(:,:,SubjInd);
    
    
    % Dirty trick in case "there was no target or distractor for a given
    % condition for a given block"
    HIT_Block_TOTAL((MISS_Block_TOTAL+HIT_Block_TOTAL)==0)=NaN;
    CORRECT_REJECTION_Block_TOTAL((FALSE_ALARM_Block_TOTAL+CORRECT_REJECTION_Block_TOTAL)==0)=NaN;
    FALSE_ALARM_Block_TOTAL((FALSE_ALARM_Block_TOTAL+CORRECT_REJECTION_Block_TOTAL)==0)=NaN;
    
    
    
    
    
    fprintf('Subject %i average', SubjInd)
    
    GroupResults(SubjInd).HitRate(:,:,:) = HIT_TOTAL./(HIT_TOTAL+MISS_TOTAL); %#ok<*SAGROW>
    GroupResults(SubjInd).FARate(:,:,:) = FALSE_ALARM_TOTAL./(FALSE_ALARM_TOTAL+CORRECT_REJECTION_TOTAL); %#ok<*SAGROW>
    GroupResults(SubjInd).CorrectRejectionRate(:,:,:) = CORRECT_REJECTION_TOTAL./(CORRECT_REJECTION_TOTAL+FALSE_ALARM_TOTAL);
    GroupResults(SubjInd).HitRateBlock(:,:,:) = HIT_Block_TOTAL./(HIT_Block_TOTAL+MISS_Block_TOTAL); %#ok<*SAGROW>
    GroupResults(SubjInd).CorrectRejectionRateBlock(:,:,:) = CORRECT_REJECTION_Block_TOTAL./(CORRECT_REJECTION_Block_TOTAL+FALSE_ALARM_Block_TOTAL);
    GroupResults(SubjInd).FARateBlock(:,:,:) = FALSE_ALARM_Block_TOTAL./(CORRECT_REJECTION_Block_TOTAL+FALSE_ALARM_Block_TOTAL);
    
    HIT_TOTAL = nansum(HIT_TOTAL,3);
    MISS_TOTAL = nansum(MISS_TOTAL,3);
    FALSE_ALARM_TOTAL = nansum(FALSE_ALARM_TOTAL,3);
    CORRECT_REJECTION_TOTAL = nansum(CORRECT_REJECTION_TOTAL,3);
    
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
            if isnan(D_prime_TOTAL(i,j))
                warning('NaN')
            end
            Accuracy_TOTAL(i,j) = (HIT_TOTAL(i,j) + CORRECT_REJECTION_TOTAL(i,j)) / (HIT_TOTAL(i,j)+MISS_TOTAL(i,j)+FALSE_ALARM_TOTAL(i,j)+CORRECT_REJECTION_TOTAL(i,j));
            
        end
    end
    
    D_prime_TOTAL;
    Accuracy_TOTAL;
    
    GroupResults(SubjInd).MISS_Block_TOTAL = MISS_Block_TOTAL;
    
    GroupResults(SubjInd).HIT_TOTAL = HIT_TOTAL;
    GroupResults(SubjInd).CORRECT_REJECTION_TOTAL = CORRECT_REJECTION_TOTAL;
    GroupResults(SubjInd).MISS_TOTAL = MISS_TOTAL;
    GroupResults(SubjInd).FALSE_ALARM_TOTAL = FALSE_ALARM_TOTAL;
    
    GroupResults(SubjInd).D_prime_TOTAL = D_prime_TOTAL;
    GroupResults(SubjInd).Accuracy_TOTAL = Accuracy_TOTAL;
    
    for i=1:3
        GroupResults(SubjInd).AttModInd.Unscaled(i,1) = sum(HIT_TOTAL(i,1)) - sum(FALSE_ALARM_TOTAL(i,2));
        GroupResults(SubjInd).AttModInd.Unscaled(i,2) = sum(HIT_TOTAL(i,2)) - sum(FALSE_ALARM_TOTAL(i,1));
        GroupResults(SubjInd).AttModInd.Unscaled(i,3) = sum(sum(HIT_TOTAL(i,1:2))) - sum(sum(FALSE_ALARM_TOTAL(i,1:2)));
        
        GroupResults(SubjInd).AttModInd.Scaled(i,1) = sum(HIT_TOTAL(i,1))/ (sum(HIT_TOTAL(i,1))+sum(MISS_TOTAL(i,1))) ...
            - sum(FALSE_ALARM_TOTAL(i,1)) / ( sum(FALSE_ALARM_TOTAL(i,1)) + sum(CORRECT_REJECTION_TOTAL(i,1)));
        GroupResults(SubjInd).AttModInd.Scaled(i,2) = sum(HIT_TOTAL(i,2))/(sum(HIT_TOTAL(i,2))+sum(MISS_TOTAL(i,2))) ...
            - sum(FALSE_ALARM_TOTAL(i,2)) / ( sum(FALSE_ALARM_TOTAL(i,2)) + sum(CORRECT_REJECTION_TOTAL(i,2)));
        GroupResults(SubjInd).AttModInd.Scaled(i,3) = sum(sum(HIT_TOTAL(i,1:2))) /(sum(sum(HIT_TOTAL(i,1:2)))+sum(sum(MISS_TOTAL(i,1:2)))) ...
            - sum(sum(FALSE_ALARM_TOTAL(i,1:2)))/(sum(sum(FALSE_ALARM_TOTAL(i,1:2))) + sum(sum(CORRECT_REJECTION_TOTAL(i,1:2))));
    end
    
    cd(StartDirectory)
    
    fprintf('\n\n')
    
end

for iSubj=1:size(GroupResults,2)
    A=isnan(GroupResults(iSubj).HitRateBlock(1:3,1:2,:));
    A=sum(A(:))/numel(A);
    B=isnan(GroupResults(iSubj).CorrectRejectionRateBlock(1:3,1:2,:));
    B=sum(B(:))/numel(B);
    C(iSubj,1:2) = [A B];
end

mkdir(fullfile(StartDirectory,'Figures','Behavioral'))
save(fullfile(StartDirectory,'Figures','Behavioral','BehavioralResults.mat'), 'GroupResults', 'SubjectList')


%% Hits and FA

close all

FontSize = 12;

figure('name', 'Hits and FA', 'position', [100, 100, 750, 500], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    Hit_tmp(:,:,iSubj) = mean(GroupResults(iSubj).HitRate,3);
    FA_tmp(:,:,iSubj) = mean(GroupResults(iSubj).FARate,3);
end

N = size(Hit_tmp,3);

Hit_tmp_mean = nanmean(Hit_tmp,3);
Hit_tmp_sem = nanstd(Hit_tmp,0,3)/N^.5;

FA_tmp_mean = nanmean(FA_tmp,3);
FA_tmp_sem = nanstd(FA_tmp,0,3)/N^.5;

for AttCdt = 1:2
    
    subplot(2, 1, AttCdt)
    hold on
    grid on
    
    for iSubj=1:size(GroupResults,2)
        plot(1:3, Hit_tmp(1:3, AttCdt, iSubj), '-o', ...
            'Color', [.5, .5, 1], 'MarkerFaceColor', [.5, .5, 1],...
            'MarkerSize', 4)
        
        plot(1:3, FA_tmp(1:3, AttCdt, iSubj), '-o', ...
            'Color', [1, .5, .5], 'MarkerFaceColor', [1, .5, .5],...
            'MarkerSize', 4)
    end
    
    errorbar([1:3]-.1, Hit_tmp_mean(1:3, AttCdt), Hit_tmp_sem(1:3, AttCdt), ...
        '-ob', 'linewidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5)
    errorbar([1:3]+.1, FA_tmp_mean(1:3, AttCdt), FA_tmp_sem(1:3, AttCdt), ...
        '-or', 'linewidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 5)
    
    set(gca, 'xtick',1:3, ...
        'xticklabel', {'Auditory','Visual', 'Audiovisual'}, ...
        'ytick',0:.2:1, ...
        'yticklabel',0:20:100, ...
        'fontsize',FontSize);
    
    axis([0.7 3.2 0 1])
    
end

subplot(2, 1, 1)
t=title('Auditory attention');
set(t,'fontsize',FontSize+2);

subplot(2, 1, 2)
t=title('Visual attention');
set(t,'fontsize',FontSize+2);


% print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'HitsAndFalseAlarm.tif'), '-dtiff')


return


%% Hits and FA

close all

Scatter=linspace(0.2,0.6,size(GroupResults,2));

FontSize = 10;

figure('name', 'Hits and FA', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    Hit_tmp(:,:,iSubj) = mean(GroupResults(iSubj).HitRate,3);
    FA_tmp(:,:,iSubj) = mean(GroupResults(iSubj).FARate,3);
end

N = size(Hit_tmp,3);

Hit_tmp_mean = nanmean(Hit_tmp,3);
Hit_tmp_sem = nanstd(Hit_tmp,0,3)/N^.5;

FA_tmp_mean = nanmean(FA_tmp,3);
FA_tmp_sem = nanstd(FA_tmp,0,3)/N^.5;

iSubplot = 1;
for iRow=1:2
    for iCol=1:3
        
        subplot(2, 3, iSubplot)
        hold on
        grid on
        axis([0 size(GroupResults,2)+.5 0 1])
        
        errorbar(0.35, Hit_tmp_mean(iCol,iRow), Hit_tmp_sem(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        errorbar(1.25, FA_tmp_mean(iCol,iRow), FA_tmp_sem(iCol,iRow), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
        
        for iSubj=1:size(GroupResults,2)
            plot(0.35+Scatter(iSubj),Hit_tmp(iCol,iRow,iSubj), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            
            plot(1.25+Scatter(iSubj),FA_tmp(iCol,iRow,iSubj), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
        end
        axis([0.2 2 0 1])
        set(gca, 'xtick',[0.35 1.25], 'xticklabel', {'Hits','False alarms'}, 'ytick', 0:.2:1, 'fontsize',FontSize);
        iSubplot = iSubplot + 1;
    end
end

subplot(2, 3, 1)
t=title('Auditory stimulation');
set(t,'fontsize',FontSize+2);
t=ylabel('Auditory attention');
set(t,'fontsize',FontSize+2);

subplot(2, 3, 4)
t=ylabel('Visual attention');
set(t,'fontsize',FontSize+2);

subplot(2, 3, 2)
t=title('Visual stimulation');
set(t,'fontsize',FontSize+2);

subplot(2, 3, 3)
t=title('AudioVisual stimulation');
set(t,'fontsize',FontSize+2);

% print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'HitsAndFalseAlarm.tif'), '-dtiff')

clear TEMP TEMP2 TEMP3

return


%% Save tables D prime
MEAN = nanmean(cat(3,GroupResults.D_prime_TOTAL),3);
SEM = nansem(cat(3,GroupResults.D_prime_TOTAL),3);

SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','d_prime.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'D prime,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,1), SEM(i,1));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,2), SEM(i,2));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Total,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,3), SEM(i,3));
end

fclose (fid);



SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','d_prime_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:numel(GroupResults)
    fprintf(fid, '%f,', GroupResults(iSubj).D_prime_TOTAL([1:3 5:7]));
    fprintf (fid, '\n');
end

fclose (fid);



%% Save tables Accuracy
MEAN = nanmean(cat(3,GroupResults.Accuracy_TOTAL),3);
SEM = nansem(cat(3,GroupResults.Accuracy_TOTAL),3);

SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','accuracy.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Accuracy,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,1), SEM(i,1));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,2), SEM(i,2));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Total,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,3), SEM(i,3));
end

fclose (fid);


SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','accuracy_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:numel(GroupResults)
    fprintf(fid, '%f,', GroupResults(iSubj).Accuracy_TOTAL([1:3 5:7]));
    fprintf (fid, '\n');
end

fclose (fid);



%% Save tables Hit Rate
MEAN = nanmean(cat(3,GroupResults.HitRate),3);
SEM = nansem(cat(3,GroupResults.HitRate),3);

SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','hit_rate.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'D prime,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,1), SEM(i,1));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,2), SEM(i,2));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Total,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,3), SEM(i,3));
end

fclose (fid);



SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','hit_rate_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:numel(GroupResults)
    fprintf(fid, '%f,', GroupResults(iSubj).HitRate([1:3 5:7]));
    fprintf (fid, '\n');
end

fclose (fid);



%% Save tables CR Rate
MEAN = nanmean(cat(3,GroupResults.CorrectRejectionRate),3);
SEM = nansem(cat(3,GroupResults.CorrectRejectionRate),3);

SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','correct_rejection_rate.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'D prime,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,1), SEM(i,1));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,2), SEM(i,2));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Total,');
for i=1:4
    fprintf (fid, '%f, (%f),', MEAN(i,3), SEM(i,3));
end

fclose (fid);



SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','correct_rejection_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:numel(GroupResults)
    fprintf(fid, '%f,', GroupResults(iSubj).CorrectRejectionRate([1:3 5:7]));
    fprintf (fid, '\n');
end

fclose (fid);



%% Save tables RT Hit
SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','RT_Hit_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:numel(GroupResults)
    TEMP = cellfun(@median, GroupResults(iSubj).RT_HIT);
    RT_Hit_all(:,:,iSubj) = TEMP;
    fprintf(fid, '%f,', TEMP(:));
    fprintf (fid, '\n');
end

fclose (fid);


MEDIAN = nanmedian(RT_Hit_all,3);
SEM = nansem(RT_Hit_all,3);

SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','RT_Hit.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'Accuracy,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:3
    fprintf (fid, '%f, (%f),', MEDIAN(i,1), SEM(i,1));
end
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=1:3
    fprintf (fid, '%f, (%f),', MEDIAN(i,2), SEM(i,2));
end

fclose (fid);


%% Targets

close all

FontSize = 6;

figure('name', 'Attention Modulation Indices', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    A_Targets(iSubj,:) = [...
        GroupResults(iSubj).HIT_TOTAL(1:3,1)'...
        GroupResults(iSubj).FALSE_ALARM_TOTAL(1:3,2)'];
    
    V_Targets(iSubj,:) = [...
        GroupResults(iSubj).FALSE_ALARM_TOTAL(1:3,1)'...
        GroupResults(iSubj).HIT_TOTAL(1:3,2)'];
end

subplot(2,2,1)
hold on ; grid on

for iStim=1:3
    
    errorbar(iStim, nanmean(A_Targets(:,iStim,:)), nansem(A_Targets(:,iStim,:)), ...
        '.b', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    
    h = plotSpread(A_Targets(:,iStim,:), 'distributionIdx', ones(size(A_Targets(:,iStim,:))), ...
        'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
        'xValues', iStim+.4, 'binWidth', 0.1, 'spreadWidth', 0.8);
    set(h{1}, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1)
end

axis([0.8 3.5 -0.5 16.5])
set(gca, 'xtick',1.1:3.1, 'xticklabel', {'Audio' 'Visual' 'Audiovisual'}, 'ytick', 0:2:100, 'fontsize',FontSize+2, 'xgrid', 'off');

t=title('Auditory targets');
set(t,'fontsize',FontSize+4);
t=ylabel('Auditory attention');
set(t,'fontsize',FontSize+4);
t=xlabel('Stimulus');
set(t,'fontsize',FontSize+2);


subplot(2,2,3)
hold on; grid on

for iStim=1:3
    
    errorbar(iStim, nanmean(A_Targets(:,iStim+3,:)), nansem(A_Targets(:,iStim+3,:)), ...
        '.b', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    
    h = plotSpread(A_Targets(:,iStim+3,:), 'distributionIdx', ones(size(A_Targets(:,iStim,:))), ...
        'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
        'xValues', iStim+.4, 'binWidth', 0.1, 'spreadWidth', 0.8);
    set(h{1}, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1)
    
end

axis([0.8 3.5 -0.5 16.5])
set(gca, 'xtick',1.1:3.1, 'xticklabel', {'Audio' 'Visual' 'Audiovisual'}, 'ytick', 0:2:100, 'fontsize',FontSize+2, 'xgrid', 'off');

t=ylabel('Visual attention');
set(t,'fontsize',FontSize+4);
t=xlabel('Stimulus');
set(t,'fontsize',FontSize+2);



subplot(2,2,2)
hold on ; grid on

for iStim=1:3
    
    errorbar(iStim, nanmean(V_Targets(:,iStim,:)), nansem(V_Targets(:,iStim,:)), ...
        '.b', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    
    h = plotSpread(V_Targets(:,iStim,:), 'distributionIdx', ones(size(V_Targets(:,iStim,:))), ...
        'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
        'xValues', iStim+.4, 'binWidth', 0.1, 'spreadWidth', 0.8);
    set(h{1}, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1)
end

axis([0.8 3.5 -0.5 16.5])
set(gca, 'xtick',1.1:3.1, 'xticklabel', {'Audio' 'Visual' 'Audiovisual'}, 'ytick', 0:2:100, 'fontsize',FontSize+2, 'xgrid', 'off');

t=title('Visual targets');
set(t,'fontsize',FontSize+4);
t=xlabel('Stimulus');
set(t,'fontsize',FontSize+2);


subplot(2,2,4)
hold on; grid on

for iStim=1:3
    
    errorbar(iStim, nanmean(V_Targets(:,iStim+3,:)), nansem(V_Targets(:,iStim+3,:)), ...
        '.b', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    
    h = plotSpread(V_Targets(:,iStim+3,:), 'distributionIdx', ones(size(V_Targets(:,iStim,:))), ...
        'distributionMarkers',{'.'},'distributionColors',{'k'}, ...
        'xValues', iStim+.4, 'binWidth', 0.1, 'spreadWidth', 0.8);
    set(h{1}, 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1)
    
end

axis([0.8 3.5 -0.5 16.5])
set(gca, 'xtick',1.1:3.1, 'xticklabel', {'Audio' 'Visual' 'Audiovisual'}, 'ytick', 0:2:100, 'fontsize',FontSize+2, 'xgrid', 'off');

t=xlabel('Stimulus');
set(t,'fontsize',FontSize+2);

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'Targets.tif'), '-dtiff')


clear TEMP TEMP2 TEMP3



% Save tables A target and V targets
SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','Targets.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'A Targets,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:3
    fprintf (fid, '%f, (%f),', nanmean(A_Targets(:,i)), nansem(A_Targets(:,i)));
end
fprintf (fid, '%f, (%f),', nanmean(sum(A_Targets(:,1:3),2)), nansem(sum(A_Targets(:,1:3)),2));
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=4:6
    fprintf (fid, '%f, (%f),', nanmean(A_Targets(:,i)), nansem(A_Targets(:,i)));
end
fprintf (fid, '%f, (%f),', nanmean(sum(A_Targets(:,4:6),2)), nansem(sum(A_Targets(:,4:6)),2));
fprintf (fid, '\n\n');


fprintf (fid, '%s', 'V Targets,','stimulus,');
fprintf (fid, '\n');
fprintf (fid, '%s', ',','Audio,',',','Visual,',',','Audio-Visual,',',','Total');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to audio,');
for i=1:3
    fprintf (fid, '%f, (%f),', nanmean(V_Targets(:,i)), nansem(V_Targets(:,i)));
end
fprintf (fid, '%f, (%f),', nanmean(sum(V_Targets(:,1:3),2)), nansem(sum(V_Targets(:,1:3)),2));
fprintf (fid, '\n');
fprintf (fid, '%s', 'Attend to visual,');
for i=4:6
    fprintf (fid, '%f, (%f),', nanmean(V_Targets(:,i)), nansem(V_Targets(:,i)));
end
fprintf (fid, '%f, (%f),', nanmean(sum(V_Targets(:,4:6),2)), nansem(sum(V_Targets(:,4:6)),2));
fprintf (fid, '\n\n');

fclose (fid);



SavedTxt = fullfile('/data','AV_Integration_2','Figures','Behavioral','Targets_subj.csv');
fid = fopen (SavedTxt, 'w');

fprintf (fid, '%s', 'A Targets,');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:size(A_Targets,1)
    fprintf(fid, '%f,', A_Targets(iSubj,:));
    fprintf (fid, '\n');
end


fprintf (fid, '\n\n');
fprintf (fid, '%s', 'V Targets,');
fprintf (fid, '\n');
fprintf (fid, '%s', 'Stim_A-Att_A,', 'Stim_V-Att_A,', 'Stim_AV-Att_A,', 'Stim_A-Att_V,', 'Stim_V-Att_V,', 'Stim_AV-Att_V,');
fprintf (fid, '\n');
for iSubj=1:size(V_Targets,1)
    fprintf(fid, '%f,', V_Targets(iSubj,:));
    fprintf (fid, '\n');
end

fclose (fid);


%% Attention Modulation Indices

close all

Scatter=linspace(0.2,0.6,size(GroupResults,2));

FontSize = 6;

figure('name', 'Attention Modulation Indices', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    AttModInd.Unscaled(:,:,iSubj) = GroupResults(iSubj).AttModInd.Unscaled;
    AttModInd.Scaled(:,:,iSubj) = GroupResults(iSubj).AttModInd.Scaled;
end

iSubplot = 1;
for iCol=1:3
    
    subplot(2,3,iSubplot)
    hold on
    grid on
    %     axis([0 size(GroupResults,2)+.5 0 1])
    
    for iStim=1:3
        
        errorbar(iStim, nanmean(AttModInd.Unscaled(iStim,iCol,:),3), nansem(AttModInd.Unscaled(iStim,iCol,:),3), ...
            'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        
        for iSubj=1:size(GroupResults,2)
            plot(iStim+Scatter(iSubj),AttModInd.Unscaled(iStim,iCol,iSubj), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
        end
        
    end
    
    axis([0.5 iStim+Scatter(end)+.5 0 32])
    set(gca, 'xtick',1:3, 'xticklabel', {'A' 'V' 'AV'}, 'ytick', 0:2:100, 'fontsize',FontSize+2);
    iSubplot = iSubplot + 1;
    
end

iSubplot = 4;
for iCol=1:3
    
    subplot(2,3,iSubplot)
    hold on
    grid on
    
    for iStim=1:3
        
        errorbar(iStim, nanmean(AttModInd.Scaled(iStim,iCol,:),3), nansem(AttModInd.Scaled(iStim,iCol,:),3), ...
            'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        
        for iSubj=1:size(GroupResults,2)
            plot(iStim+Scatter(iSubj),AttModInd.Scaled(iStim,iCol,iSubj), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
        end
    end
    
    axis([0.5 iStim+Scatter(end)+.5 0 1])
    set(gca, 'xtick',1:3, 'xticklabel', {'A' 'V' 'AV'}, 'ytick', 0:.1:1, 'fontsize',FontSize+2);
    iSubplot = iSubplot + 1;
    
end



subplot(2,3,1)
t=title('Unscaled: audio att');
set(t,'fontsize',FontSize+2);
subplot(2,3,2)
t=title('Unscaled: visual att');
set(t,'fontsize',FontSize+2);
subplot(2,3,3)
t=title('Unscaled: all');
set(t,'fontsize',FontSize+2);

subplot(2,3,4)
t=title('Scaled: audio att');
set(t,'fontsize',FontSize+2);
subplot(2,3,5)
t=title('Scaled: visual att');
set(t,'fontsize',FontSize+2);
subplot(2,3,6)
t=title('Scaled: all');
set(t,'fontsize',FontSize+2);

mtit('Attention Modulation Indices', 'xoff',0,'yoff',.025)

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'AttentionModulationIndices.tif'), '-dtiff')


clear TEMP TEMP2 TEMP3


%% False alarms

close all

TEMP = [];

FontSize = 6;

figure('name', 'FalseAlarms', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = mean(GroupResults(iSubj).FARate,3);
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            hold on
            grid on
            if length(GroupResults(iSubj).FARate(iCol,iRow,:))~=4
                SubjectList(iSubj)
                iRow
                iCol
            end
            plot(repmat(iSubj, length(GroupResults(iSubj).FARate(iCol,iRow,:)), 1),...
                squeeze(GroupResults(iSubj).FARate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = nanmean(TEMP,3);
TEMP3 = nanstd(TEMP,3);

iSubplot = 1;
for iRow=1:3
    for iCol=1:4
        
        subplot(3,4,iSubplot)
        axis([0 size(GroupResults,2)+.5 0 1])
        errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:1, 'fontsize',FontSize-2);
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

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'FalseAlarms.tif'), '-dtiff')

clear TEMP TEMP2 TEMP3


%% Accuracy

close all

TEMP = [];

FontSize = 6;

figure('name', 'Accuracy', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = mean(GroupResults(iSubj).Accuracy,3);
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            hold on
            grid on
            if length(GroupResults(iSubj).Accuracy(iCol,iRow,:))~=4
                SubjectList(iSubj)
                iRow
                iCol
            end
            plot(repmat(iSubj, length(GroupResults(iSubj).Accuracy(iCol,iRow,:)), 1),...
                squeeze(GroupResults(iSubj).Accuracy(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = nanmean(TEMP,3);
TEMP3 = nanstd(TEMP,3);

iSubplot = 1;
for iRow=1:3
    for iCol=1:4
        
        subplot(3,4,iSubplot)
        axis([0 size(GroupResults,2)+.5 0 1])
        errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', 0:.2:1, 'fontsize',FontSize-2);
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

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral', 'Accuracy.tif'), '-dtiff')

clear TEMP TEMP2 TEMP3


%% D prime

TEMP = [];

figure('name', 'D prime', 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])

for iSubj=1:size(GroupResults,2)
    TEMP(:,:,iSubj) = mean(GroupResults(iSubj).d,3);
    MIN(iSubj) = min(GroupResults(iSubj).d(:));
    MAX(iSubj) = max(GroupResults(iSubj).d(:));
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            hold on
            grid on
            if length(GroupResults(iSubj).d(iCol,iRow,:))~=4
                SubjectList(iSubj)
                iRow
                iCol
            end
            plot(repmat(iSubj, length(GroupResults(iSubj).d(iCol,iRow,:)), 1), ...
                squeeze(GroupResults(iSubj).d(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = nanmean(TEMP,3);
TEMP3 = nanstd(TEMP,3);

MIN = min(MIN)-.1;
MAX = max(MAX)+.1;

iSubplot = 1;
for iRow=1:3
    for iCol=1:4
        
        subplot(3,4,iSubplot)
        plot([0 size(GroupResults,2)+.5],[0 0], '--k')
        axis([0 size(GroupResults,2)+.5 MIN MAX])
        errorbar(0.35, TEMP2(iCol,iRow), TEMP3(iCol,iRow), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
        set(gca, 'xtick',1:size(GroupResults,2), 'xticklabel', SubjectList, 'ytick', -1:1:3,'fontsize',FontSize-2);
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

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',  'D_Prime.tif'), '-dtiff')


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


print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',  'ReactionTime.tif'), '-dtiff')


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
            if length(GroupResults(iSubj).HitRate(iCol,iRow,:))~=4
                SubjectList(iSubj)
                iRow
                iCol
            end
            plot(repmat(iSubj, length(GroupResults(iSubj).HitRate(iCol,iRow,:)), 1), ...
                squeeze(GroupResults(iSubj).HitRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = nanmean(TEMP,3);
TEMP3 = nanstd(TEMP,3);

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

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',  'HitRate.tif'), '-dtiff')


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
            squeeze(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:))'
            if length(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:))~=4
                SubjectList(iSubj)
                iRow
                iCol
                
            end
            plot(repmat(iSubj, length(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), 1), ...
                squeeze(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            iSubplot = iSubplot + 1;
        end
    end
    
end

TEMP2 = nanmean(TEMP,3);
TEMP3 = nanstd(TEMP,3);

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

print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',  'CorrectRejectionRate.tif'), '-dtiff')


%% CR and hit rate per session/subject

close all

FontSize=10;

for iSubj=1:size(GroupResults,2)
    
    figure('name', ['Hit and CR Rate for subject ' SubjectList(iSubj,:)], ...
        'position',  [100, 100, 1500, 1000], 'Color', [1 1 1], ...
        'visible', 'off')
    
    iSubplot = 1;
    for iRow=1:3
        for iCol=1:4
            
            subplot(3,4,iSubplot)
            
            hold on
            grid on
            
            
            plot(1:size(GroupResults(iSubj).HitRate,3), squeeze(GroupResults(iSubj).HitRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            
            plot(.1+[1:size(GroupResults(iSubj).HitRate,3)], squeeze(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), ...
                'MarkerSize', 4)
            
            errorbar(0.35, mean(GroupResults(iSubj).HitRate(iCol,iRow,:)),...
                std(GroupResults(iSubj).HitRate(iCol,iRow,:)), ...
                'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
            
            errorbar(0.45, mean(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)),...
                std(GroupResults(iSubj).CorrectRejectionRate(iCol,iRow,:)), ...
                'ob', 'MarkerSize', 5)
            
            if iSubplot==12
                legend({'Hit rate';'CR rate'; 'mean hit rate'; 'mean CR rate'}', ...
                    'location', 'SouthEast', 'fontsize', FontSize-6)
            end
            
            plot([0 size(GroupResults,2)+.5],[0.7 0.7], '--k')
            
            axis([0 size(GroupResults(iSubj).HitRate,3)+.5 0 1])
            
            set(gca, 'xtick',1:size(GroupResults(iSubj).HitRate,3), ...
                'xticklabel', 1:size(GroupResults(iSubj).HitRate,3), ...
                'ytick', 0:.2:1,'fontsize',FontSize);
            
            t=xlabel('Session');
            set(t,'fontsize',FontSize);
            
            iSubplot = iSubplot + 1;
        end
    end
    
    
    subplot(3,4,1)
    t=title('Auditory stimulation');
    set(t,'fontsize',FontSize-2);
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
    set(t,'fontsize',FontSize-2);
    
    subplot(3,4,3)
    t=title('AudioVisual stimulation');
    set(t,'fontsize',FontSize-2);
    
    subplot(3,4,4)
    t=title('All');
    set(t,'fontsize',FontSize-2);
    
    print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',...
        ['Hit&CR_RateSubject_' SubjectList(iSubj,:) '.tif']), '-dtiff')
    
    
end


%% CR and hit rate per session/subject

close all

FontSize=10;

for iSubj=1:size(GroupResults,2)
    
    figure('name', ['Hit and CR Rate for subject ' SubjectList(iSubj,:)], 'position',  [100, 100, 1500, 1000], 'Color', [1 1 1])
    
    iSubplot = 1;
    for iRow=1:2
        for iCol=1:3
            
            subplot(2,3,iSubplot)
            
            hold on
            grid on
            
            
            plot(1:size(GroupResults(iSubj).HitRateBlock,3), ...
                squeeze(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
                'MarkerSize', 4)
            
            plot(.2+[1:size(GroupResults(iSubj).HitRateBlock,3)], ...
                squeeze(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)), ' o', ...
                'Color', COLOR_Subject(iSubj,:), ...
                'MarkerSize', 4)
            
            errorbar(0.35, mean(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)),...
                std(GroupResults(iSubj).HitRateBlock(iCol,iRow,:)), ...
                'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
            
            errorbar(0.45, mean(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)),...
                std(GroupResults(iSubj).CorrectRejectionRateBlock(iCol,iRow,:)), ...
                'ob', 'MarkerSize', 5)
            
            if iSubplot==6
                legend({'Hit rate';'CR rate'; 'mean hit rate'; 'mean CR rate'}', ...
                    'location', 'SouthEast', 'fontsize', FontSize-6)
            end
            
            plot([0 size(GroupResults(iSubj).HitRateBlock,3)+.5],[0.7 0.7], '--k')
            plot([3.5 3.5],[0 1], '-k', 'linewidth', 2)
            plot([6.5 6.5],[0 1], '-k', 'linewidth', 2)
            plot([9.5 9.5],[0 1], '-k', 'linewidth', 2)
            
            axis([0 size(GroupResults(iSubj).HitRateBlock,3)+.5 0 1])
            
            set(gca, 'xtick',1:size(GroupResults(iSubj).HitRateBlock,3), ...
                'xticklabel', 1:size(GroupResults(iSubj).HitRateBlock,3), ...
                'ytick', 0:.2:1,'fontsize',FontSize);
            
            t=xlabel('Blocks');
            set(t,'fontsize',FontSize);
            
            iSubplot = iSubplot + 1;
        end
    end
    
    
    subplot(2,3,1)
    t=title('Auditory stimulation');
    set(t,'fontsize',FontSize-2);
    t=ylabel('Auditory attention');
    set(t,'fontsize',FontSize);
    
    subplot(2,3,4)
    t=ylabel('Visual attention');
    set(t,'fontsize',FontSize);
    
    subplot(2,3,2)
    t=title('Visual stimulation');
    set(t,'fontsize',FontSize-2);
    
    subplot(2,3,3)
    t=title('AudioVisual stimulation');
    set(t,'fontsize',FontSize-2);
    
    print(gcf, fullfile('/data','AV_Integration_2','Figures','Behavioral',...
        ['Hit&CR_Rate_PerBlock_Subject_' SubjectList(iSubj,:) '.tif']), '-dtiff')
    
    
end

