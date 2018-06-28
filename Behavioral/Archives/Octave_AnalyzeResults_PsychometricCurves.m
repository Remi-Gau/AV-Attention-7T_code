clc; clear all; close all

addpath(fullfile(pwd, 'Palamedes'));

Subject2Check = 16;

% message = 'Parametric Bootstrap (1) or Non-Parametric Bootstrap? (2): ';
ParOrNonPar = 1; %input(message);
%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = 0:.05:1;
searchGrid.beta = logspace(1,3,100);
searchGrid.gamma = 0;  %scalar here (since fixed) but may be vector
searchGrid.lambda = 0:.05:1;  %ditto

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 1];  %1: free parameter, 0: fixed parameter

%Fit a Logistic function
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant

%Optional:
options = PAL_minimize('options');   %type PAL_minimize('options','help') for help
options.TolFun = 1e-09;     %increase required precision on LL
options.MaxIter = 10;
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 1];        %limit range for lambda
                            %(will be ignored here since lambda is not a
                            %free parameter)


IndStart = [5 ; 12];% first row of data points in txt file


StartFolder = pwd; 

SubjectList = dir('Subject*');
%IsDir = [SubjectList(:).isdir]
%SubjectList(IsDir == 0) = []
%SubjectList(1:2) = []


for SubjInd=Subject2Check

        fprintf('\n\n')
        message = sprintf('Subject %i', SubjInd);
        disp(message);

        cd (fullfile (pwd, ['Subject_' num2str(Subject2Check)], 'PsychometricCurves'))

        LogFileList = dir('Logfile*.txt');
        
        for FileInd = 1:size(LogFileList,1)
            
            if LogFileList(FileInd).name(end-22:end-20)=='111' %#ok<STCMP>
                Audio_OR_Visual = 1;
                fprintf('\n')
                message = sprintf('Auditory targets');
                disp(message);
            else
                Audio_OR_Visual = 2;
                fprintf('\n')
                message = sprintf('Visual targets');
                disp(message);                
            end
            
            % Reads the log file with the target informations
            TargetLogfile = fopen(fullfile (pwd, strcat('Onset_Events', LogFileList(FileInd).name(8:end))));
            
            % with octave
            FileContent = dlmread(TargetLogfile, '\t');
            
            fclose(TargetLogfile);
            
            
            
            if Audio_OR_Visual == 1;
                Target_Info = FileContent(:,[1 6]);
            else
                Target_Info = FileContent(:,[1 4]);
            end
            
            
            % Octave
            TEMP = [find(Target_Info(:,1)==0) ; find(Target_Info(:,1)==1) ; find(Target_Info(:,1)==2)] ;        
                     
            Target_Info(TEMP,:) = [];
            
            clear TEMP
            
            
            if Audio_OR_Visual == 1;
                LastTrial = find(Target_Info(:,1)==10, 1, 'last' );
                IndexTargets = find(Target_Info(:,1)==10)-1;
            else
                LastTrial = find(Target_Info(:,1)==11, 1, 'last' );
                IndexTargets = find(Target_Info(:,1)==11)-1;
            end
            
            IndexDisctrators = setxor(find(Target_Info(:,1)>=100), IndexTargets); 
            
            TEMP = find(Target_Info(:,1)<100);
            Target_Info(TEMP-1,2) = Target_Info(TEMP,2);
            Target_Info(TEMP,:) = [];
            
            clear TEMP
            
            
            % Reads the general logfile
            disp(LogFileList(FileInd).name)
            GeneralLogfile = fopen(fullfile (pwd, LogFileList(FileInd).name));
            FileContent = textscan(GeneralLogfile, '%s', 'delimiter','\t', 'headerlines', IndStart(1), 'returnOnError',0);
            fclose(GeneralLogfile);
            
            EOF = find(strcmp('Final_Fixation', FileContent{1,1}));
            Stim_Time = FileContent{1,1}(1:EOF(1)-3);

            clear EOF FileContent
            
            
            ToKeep = [];
            
            TEMP = find(strcmp('AudioOnly_Trial', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP
            
            TEMP = find(strcmp('VisualOnly_Trial', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP
            
            TEMP = find(strcmp('AudioVisual_Trial', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP
            
            TEMP = find(strcmp('Response', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP
            
            TEMP = find(strcmp('Visual_Target', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP
            
            TEMP = find(strcmp('Auditory_Target', Stim_Time));
            ToKeep = [ToKeep ; TEMP(:)];
            clear TEMP    
            
            Stim_Time = {Stim_Time{sort(ToKeep(:))}}';
            
            TEMP = nan(length(Stim_Time),1);
            
            TEMP(find(strcmp('AudioOnly_Trial', Stim_Time)),:) = 100;
            TEMP(find(strcmp('VisualOnly_Trial', Stim_Time)),:) = 200;
            TEMP(find(strcmp('AudioVisual_Trial', Stim_Time)),:) = 300;
            TEMP(find(strcmp('1', Stim_Time)),:) = 1;
            TEMP(find(strcmp('Response', Stim_Time)),:) = 1;
            TEMP(find(strcmp('Visual_Target', Stim_Time)),:) = 11;
            TEMP(find(strcmp('Auditory_Target', Stim_Time)),:) = 10;
            
            
            Stim_Time = TEMP;
            
            clear TEMP
            
            Stim_Time = [Stim_Time Stim_Time==1];
            Stim_Time(find(Stim_Time(:,2)==1)-1,2)=1;
            
            if Audio_OR_Visual == 1;
                Stim_Time(find(all(Stim_Time == repmat([10 1],length(Stim_Time),1),2 ) )-1 , 2 ) = 1;
                Stim_Time( all( Stim_Time == repmat([10 1],length(Stim_Time),1), 2) , : ) = [];
                Stim_Time( all( Stim_Time == repmat([1 1],length(Stim_Time),1), 2) , : ) = [];
                Stim_Time( all( Stim_Time == repmat([10 0],length(Stim_Time),1), 2) , : ) = [];
            else
                Stim_Time(find(all(Stim_Time == repmat([11 1],length(Stim_Time),1),2 ) )-1 , 2 ) = 1;
                Stim_Time( all( Stim_Time == repmat([11 1],length(Stim_Time),1), 2) , : ) = [];
                Stim_Time( all( Stim_Time == repmat([1 1],length(Stim_Time),1), 2) , : ) = [];
                Stim_Time( all( Stim_Time == repmat([11 0],length(Stim_Time),1), 2) , : ) = [];
            end


            if Audio_OR_Visual == 1;
                figure('name', strcat('Subject ', num2str(SubjInd), ' MLE Psychometric Function Fitting : Auditory targets'));
            else
                figure('name', strcat('Subject ', num2str(SubjInd), ' MLE Psychometric Function Fitting : Visual targets'));
            end
            
            
            for i = 1:4
                
                NumberCorrectPerLevel=[];
                
                switch i
                    % Gets Values for all type of trials
                    case 1
                        %fprintf('\n All trials\n')
                        LevelsUsed = unique(Target_Info(isnan(Target_Info(:,2))==0,2));
                        TrialPerLevel = hist(Target_Info(isnan(Target_Info(:,2))==0,2), LevelsUsed);

                        for LevelInd = 1:length(LevelsUsed)
                            NumberCorrectPerLevel(LevelInd) = sum(Stim_Time(Target_Info(:,2)==LevelsUsed(LevelInd),2));  %#ok<AGROW>
                        end

                    case 2
                        %fprintf('\n Auditory trials\n')
                        LevelsUsed = unique(Target_Info(all([Target_Info(:,1)==100 isnan(Target_Info(:,2))==0],2),2));
                        TrialPerLevel = hist(Target_Info( all([Target_Info(:,1)==100 isnan(Target_Info(:,2))==0],2) ,2), LevelsUsed);
                        
                        for LevelInd = 1:length(LevelsUsed) 
                            NumberCorrectPerLevel(LevelInd) = sum(Stim_Time(all([Target_Info(:,1)==100 Target_Info(:,2)==LevelsUsed(LevelInd)],2),2));  %#ok<AGROW>
                        end
                        
                    case 3
                        %fprintf('\n Visual trials\n')
                        LevelsUsed = unique(Target_Info(all([Target_Info(:,1)==200 isnan(Target_Info(:,2))==0],2),2));
                        TrialPerLevel = hist(Target_Info( all([Target_Info(:,1)==200 isnan(Target_Info(:,2))==0],2) ,2), LevelsUsed);
                        
                        for LevelInd = 1:length(LevelsUsed) 
                            NumberCorrectPerLevel(LevelInd) = sum(Stim_Time(all([Target_Info(:,1)==200 Target_Info(:,2)==LevelsUsed(LevelInd)],2),2));  %#ok<AGROW>
                        end
                        
                    case 4
                        %fprintf('\n Audio-visual trials\n')
                        LevelsUsed = unique(Target_Info(all([Target_Info(:,1)==300 isnan(Target_Info(:,2))==0],2),2));
                        TrialPerLevel = hist(Target_Info( all([Target_Info(:,1)==300 isnan(Target_Info(:,2))==0],2) ,2), LevelsUsed);
                        
                        for LevelInd = 1:length(LevelsUsed) 
                            NumberCorrectPerLevel(LevelInd) = sum(Stim_Time(all([Target_Info(:,1)==300 Target_Info(:,2)==LevelsUsed(LevelInd)],2),2));  %#ok<AGROW>
                        end                
                end

                % Psychometric curve
                LevelsUsedFineGrain=[min(LevelsUsed):max(LevelsUsed)./200:max(LevelsUsed)];

                [paramsValues LL exitflag output] = PAL_PFML_Fit(LevelsUsed',NumberCorrectPerLevel, ...
                    TrialPerLevel,searchGrid,paramsFree,PF,'searchOptions',options, ...
                    'lapseLimits',lapseLimits);
                
                %message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
                %disp(message);

                ProportionCorrectModel = PF(paramsValues, LevelsUsedFineGrain);

                subplot(4,1,i)

                plot(LevelsUsed,NumberCorrectPerLevel./TrialPerLevel,'k.','markersize',10);
                %set(gca, 'fontsize',16);
                set(gca, 'Xtick',LevelsUsed);
                axis([min(LevelsUsed) max(LevelsUsed) 0 1]);
                hold on;
                
                plot(LevelsUsedFineGrain,ProportionCorrectModel,'g-','linewidth',4);
                
                if i==4
                    if Audio_OR_Visual == 1;
                        xlabel('Attenuation index');
                    else
                        xlabel('Stimulus size');
                    end          
                end
                
            
            end
            
            subplot(4,1,1)
            if Audio_OR_Visual == 1;
                title('Auditory targets : proportion correct = f(Target intensity)');
            else
                title('Visual targets : proportion correct = f(Target size)');
            end
            ylabel('All trials');
            
            subplot(4,1,2)
            ylabel('Auditory trials');
            
            subplot(4,1,3)
            ylabel('Visual trials');
            
            subplot(4,1,4)
            ylabel('AV trials');
            

        end

        cd (StartFolder)

end
