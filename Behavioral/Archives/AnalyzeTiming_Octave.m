clc
clear all
close all

NbSubject = 2;

StartDirectory = pwd;

for SubjInd = min(NbSubject) : max(NbSubject)

    cd(strcat('Subject_', num2str(SubjInd)))
    
    cd('Runs')

    %%
    LogFileList = dir('Logfile*.txt');
    
    SOT_Events = cell(3,2,size(LogFileList,1));
    SOT_Blocks = cell(3,2,size(LogFileList,1));
    SOT_Targets = cell(1,size(LogFileList,1));
    SOT_Responses = cell(1,size(LogFileList,1));

    IndStart = 5; % first row of data points in txt file

    %%
    for FileInd = 1:size(LogFileList,1)
    
        IndexStartAudioAttention = [] ;
        IndexStartVisualAttention = [];
        IndexStartExperiment = [];
        AudioOnly_Index = cell(1,2);
        VisualOnly_Index = cell(1,2);
        AudioVisual_Index = cell(1,2);
        Visual_Target_Index = [];
        Auditory_Target_Index = [];
        Responses_Index = [];

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


        
        %%
        % Identify the two attenttional component of the run
        IndexStartAudioAttention = find(strcmp('Attend2Audio_Fixation', Stim_Time))+1;
        if isempty(IndexStartAudioAttention)
            IndexStartAudioAttention = inf;
        end

        IndexStartVisualAttention = find(strcmp('Attend2Visual_Fixation', Stim_Time))+1;
        if isempty(IndexStartVisualAttention)
            IndexStartVisualAttention = inf;
        end
        
        
        % Identify the first condition of the run
        if IndexStartAudioAttention < IndexStartVisualAttention
            Conditions = [1 2];
            IndexConditions = [IndexStartAudioAttention IndexStartVisualAttention];
        else
            Conditions = [2 1];
            IndexConditions = [IndexStartVisualAttention IndexStartAudioAttention];
        end
        
        IndexConditions;
        
        % Start XP
        IndexStartExperiment = find(strcmp('Pulse', Stim_Time))+2;
        IndexStartExperiment = sort(IndexStartExperiment);
        IndexStartExperiment = IndexStartExperiment(1);
        
        % Trials
        TEMP = find(strcmp('AudioOnly_Trial', Stim_Time))+1 ;
        AudioOnly_Index{1,Conditions(1)} = TEMP(TEMP>IndexConditions(1) & TEMP<IndexConditions(2));
        AudioOnly_Index{1,Conditions(2)} = TEMP(TEMP>IndexConditions(2));
        clear TEMP

        TEMP = find(strcmp('VisualOnly_Trial', Stim_Time))+1;        
        VisualOnly_Index{1,Conditions(1)} = TEMP(TEMP>IndexConditions(1) & TEMP<IndexConditions(2));
        VisualOnly_Index{1,Conditions(2)} = TEMP(TEMP>IndexConditions(2));
        clear TEMP

        TEMP = find(strcmp('AudioVisual_Trial', Stim_Time))+1;
        AudioVisual_Index{1,Conditions(1)} = TEMP(TEMP>IndexConditions(1) & TEMP<IndexConditions(2));
        AudioVisual_Index{1,Conditions(2)} = TEMP(TEMP>IndexConditions(2));
        clear TEMP
        
        % Targets
        Visual_Target_Index = find(strcmp('Visual_Target', Stim_Time))+1;

        Auditory_Target_Index = find(strcmp('Auditory_Target', Stim_Time))+1;
        
        % Responses
        Responses_Index = find(strcmp('Response', Stim_Time))+2;
    
    
        %  Get the actual SOTs
        SOT_Targets(1,FileInd) = (str2double(char(Stim_Time(sort([Visual_Target_Index ; Auditory_Target_Index])))) ...
                                -str2double(char(Stim_Time(IndexStartExperiment))))/10000;
                                
        SOT_Responses(1,FileInd) = (str2double(char(Stim_Time(sort([Responses_Index])))) ...
                                -str2double(char(Stim_Time(IndexStartExperiment))))/10000;                                
       
        for i = 1:2
                SOT_Events{1,i,FileInd} = (str2double(char(Stim_Time(sort(AudioOnly_Index{1,i})))) ...
                                -str2double(char(Stim_Time(IndexStartExperiment))))/10000;
                SOT_Blocks{1,i,FileInd} = SOT_Events{1,i,FileInd}(1+13*(0:length(SOT_Events{1,i,FileInd})/13-1));
                
                SOT_Events{2,i,FileInd} = (str2double(char(Stim_Time(sort(VisualOnly_Index{1,i})))) ...
                                -str2double(char(Stim_Time(IndexStartExperiment))))/10000;
                SOT_Blocks{2,i,FileInd} = SOT_Events{2,i,FileInd}(1+13*(0:length(SOT_Events{2,i,FileInd})/13-1));                                
                
                SOT_Events{3,i,FileInd} = (str2double(char(Stim_Time(sort(AudioVisual_Index{1,i})))) ...
                                -str2double(char(Stim_Time(IndexStartExperiment))))/10000;
                SOT_Blocks{3,i,FileInd} = SOT_Events{3,i,FileInd}(1+13*(0:length(SOT_Events{3,i,FileInd})/13-1));                                
        end

    end
    
    
    cd(StartDirectory)

end
