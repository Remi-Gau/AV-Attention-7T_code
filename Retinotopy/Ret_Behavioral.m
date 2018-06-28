%%
clear; clc

IndStart = 5;

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

%  Root folder definition
StartFolder = pwd;

Resp = nan(size(SubjectList,1),2); %#ok<*SAGROW>
Missed = nan(size(SubjectList,1),2);
ExtraResp = nan(size(SubjectList,1),2);
NbTrials = nan(size(SubjectList,1),2);
RT = cell(size(SubjectList,1),2);


%% Define folders, number of runs, scans per run...
for SubjInd = 1:size(SubjectList,1)
    
    %% Subject's Identity and folders
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    ResponseFolder = fullfile(SubjectFolder, 'Retinotopy', 'Behavioral');
    
    cd(ResponseFolder)
    
    LogFileList = dir('Logfile*.txt');
       
    for iFile = 1:size(LogFileList,1)
        
        cd(ResponseFolder)
        
        disp(LogFileList(iFile).name)
        
        Direction(iFile) = LogFileList(iFile).name(25);
        
        fid = fopen(fullfile (pwd, LogFileList(iFile).name));
        FileContent = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s', 'headerlines', IndStart, 'returnOnError',0);
        fclose(fid);
        
        Stim_Time{1,1} = FileContent{1,3};
        Stim_Time{1,2} = char(FileContent{1,4});
        
        % Number of volumes to discard and total number of volumes to keep
        Remove = [];
        Remove = [Remove ; find(strcmp(Stim_Time{1,1}, '30'))]; %#ok<*AGROW>
        Remove = [Remove ; find(strcmp(Stim_Time{1,1}, 'NewCycle'))];
        Remove = [Remove ; find(strcmp(Stim_Time{1,1}, 'START'))];
        Remove = [Remove ; find(strcmp(Stim_Time{1,1}, 'Instruction'))];
        Remove = [Remove ; find(strcmp(Stim_Time{1,1}, 'Final_Fixation'))];
        Remove = unique(Remove);
        
        Stim_Time{1,1}(Remove) = [];
        Stim_Time{1,2}(Remove,:) = [];
        
        Resp(SubjInd,iFile) = 0; %#ok<*SAGROW>
        Missed(SubjInd,iFile) = 0;
        ExtraResp(SubjInd,iFile)= 0;
        RT{SubjInd,iFile} = [];
        IsTrial = [];
        
        NbTrials(SubjInd,iFile) = sum(strcmp(Stim_Time{1,1}, 'Visual_Target'));
        
        for i=1:numel(Stim_Time{1,1})
            
            if strcmp(Stim_Time{1,1}(i), 'Visual_Target')
                if IsTrial
                    Missed(SubjInd,iFile) = Missed(SubjInd,iFile) + 1;
                end
                if i==numel(Stim_Time{1,1})
                    Missed(SubjInd,iFile) = Missed(SubjInd,iFile) + 1;
                end
                
                IsTrial= 1;
                StimTime = str2double(Stim_Time{1,2}(i,:));
            end
            
            if strcmp(Stim_Time{1,1}(i), '1')
                if IsTrial
                    Resp(SubjInd,iFile) = Resp(SubjInd,iFile) + 1;
                    RT{SubjInd,iFile}(end+1) = (str2double(Stim_Time{1,2}(i,:))-StimTime)/10000;
                else
                    ExtraResp(SubjInd,iFile) = ExtraResp(SubjInd,iFile) + 1; %#ok<*UNRCH>
                end
                IsTrial = 0;
            end
            
        end
        
        if NbTrials(SubjInd,iFile)~=(Resp(SubjInd,iFile)+Missed(SubjInd,iFile))
            error('We are missing some trials.')
        end
        
    end
    
end


%% Accuracy
close all

Acc = Resp./NbTrials;

FontSize = 8;

figure('name', 'Accuracy', 'position', [100, 100, 1500, 1000], 'Color', [1 1 1])

hold on
grid on

for iSubj=1:size(Acc,1)

    plot(repmat(iSubj,1,sum(~isnan(Acc(iSubj,:)))), Acc(iSubj,:) , ' o', ...
        'Color', COLOR_Subject(iSubj,:), 'MarkerFaceColor', COLOR_Subject(iSubj,:),...
        'MarkerSize', 4)

end

errorbar(0.35, mean(nanmean(Acc,2)), nanstd(nanmean(Acc,2)), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5)

axis([0 size(Acc,1)+.5 0 1])
set(gca, 'xtick',1:size(Acc,1), 'xticklabel', SubjectList, 'ytick', 0:.2:1, 'fontsize',FontSize);

t=title('Retinotopy in scanner behavioral results');
set(t,'fontsize',FontSize);
t=ylabel('Accuracy');
set(t,'fontsize',FontSize);
t=xlabel('Subjects');
set(t,'fontsize',FontSize);
