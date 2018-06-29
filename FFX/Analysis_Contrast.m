clc; clear;

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


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

% StartFolder = fullfile(pwd, '..','..');

StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

Conditions_Names = {...
    'A Stim - Auditory Attention', ...
    'V Stim - Auditory Attention', ...
    'AV Stim - Auditory Attention', ...
    'A Stim - Visual Attention', ...
    'V Stim - Visual Attention', ...
    'AV Stim - Visual Attention', ...
    %     'A Targets - Auditory Attention', ...
    %     'A Targets - Visual Attention', ...
    %     'V Targets - Auditory Attention', ...
    %     'V Targets - Visual Attention', ...
    %     'Responses - Auditory Attention', ...
    %     'Responses - Visual Attention',...
    };


for SubjInd = 1:size(SubjectList,1)
    
    ContrastsNames={};
    Contrasts = [];
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    %         AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block');
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block_Smooth');
    %     AnalysisFolder = fullfile(SubjectFolder, 'FFX_NoSmooth');
    
    
    cd(AnalysisFolder)
    
    load SPM.mat
    
    tmp = char(SPM.xX.name');
    tmp(:,1:6)=[];
    B=cellstr(tmp); clear tmp
    
    %     % SIMPLE CONTRASTS
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     for CondInd=1:length(Conditions_Names)
    %
    %         Contrasts = [Contrasts ; strcmp([Conditions_Names{CondInd} '*bf(1)'], B)']; %#ok<*AGROW>
    %         ContrastsNames{end+1}=strcat(Conditions_Names{CondInd}, ' > Baseline');  %#ok<*SAGROW>
    %
    %         Contrasts = [Contrasts ; -1*strcmp([Conditions_Names{CondInd} '*bf(1)'], B)'];
    %         ContrastsNames{end+1}=strcat(Conditions_Names{CondInd}, ' < Baseline');
    %         TEMP(strcmp([Conditions_Names{CondInd} '*bf(1)'], B)) = 1;
    %     end
    
    AStim = ~cellfun('isempty',strfind(B,'A Stim'));
    VStim = ~cellfun('isempty',strfind(B,'V Stim'));
    AVStim = ~cellfun('isempty',strfind(B,'AV Stim'));
    
    TimeDer = ~cellfun('isempty',strfind(B,'*bf(2)'));
    
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AStim) = 1;
    TEMP(TimeDer)=0;
    Contrasts = [Contrasts ; TEMP];
    ContrastsNames{end+1}='A_Stim > Baseline';  %#ok<*SAGROW>
    Contrasts = [Contrasts ; -1*TEMP];
    ContrastsNames{end+1}='A_Stim < Baseline';
    
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(AVStim) = 1;
    %     TEMP(TimeDer)=0;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='AV Stim > Baseline';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='AV Stim < Baseline';
    
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(VStim) = 1;
    %     TEMP(AVStim) = 0;
    TEMP(TimeDer)=0;
    Contrasts = [Contrasts ; TEMP];
    ContrastsNames{end+1}='V_AV_Stim > Baseline';
    Contrasts = [Contrasts ; -1*TEMP];
    ContrastsNames{end+1}='V_AV_Stim < Baseline';
    
    
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(all([~AVStim VStim],2)) = 1;
    %     TEMP(AVStim) = 0;
    TEMP(TimeDer)=0;
    Contrasts = [Contrasts ; TEMP];
    ContrastsNames{end+1}='V_Stim > Baseline';
    Contrasts = [Contrasts ; -1*TEMP];
    ContrastsNames{end+1}='V_Stim < Baseline';
    
    
    
    %     % ALL
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}= 'All > Baseline';
    %
    %     Contrasts = [Contrasts ; -TEMP];
    %     ContrastsNames{end+1}= 'All < Baseline';
    %
    %
    %     % DIFFERENTIAL CONTRASTS
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(strcmp('A Stim - Auditory Attention*bf(1)', B)) = 1;
    %     TEMP(strcmp('A Stim - Visual Attention*bf(1)', B)) = -1;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='A Stim - A Att > A Stim - V Att';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='A Stim - V Att > A Stim - A Att';
    %
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(strcmp('V Stim - Auditory Attention*bf(1)', B)) = 1;
    %     TEMP(strcmp('V Stim - Visual Attention*bf(1)', B)) = -1;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='V Stim - A Att > V Stim - V Att';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='V Stim - V Att > V Stim - A Att';
    %
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(strcmp('AV Stim - Auditory Attention*bf(1)', B)) = 1;
    %     TEMP(strcmp('AV Stim - Visual Attention*bf(1)', B)) = -1;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='AV Stim - A Att > AV Stim - V Att';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='AV Stim - V Att > AV Stim - A Att';
    %
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(strcmp('AV Stim - Auditory Attention*bf(1)', B)) = 1;
    %     TEMP(strcmp('A Stim - Auditory Attention*bf(1)', B)) = -1;
    %     TEMP(strcmp('V Stim - Auditory Attention*bf(1)', B)) = -1;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='(AV > A + V) - A Att';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='(AV < A + V) - A Att';
    %
    %     TEMP = zeros(1, size(SPM.xX.X,2));
    %     TEMP(strcmp('AV Stim - Visual Attention*bf(1)', B)) = 1;
    %     TEMP(strcmp('A Stim - Visual Attention*bf(1)', B)) = -1;
    %     TEMP(strcmp('V Stim - Visual Attention*bf(1)', B)) = -1;
    %     Contrasts = [Contrasts ; TEMP];
    %     ContrastsNames{end+1}='(AV > A + V) - V Att';
    %     Contrasts = [Contrasts ; -1*TEMP];
    %     ContrastsNames{end+1}='(AV < A + V) - V Att';
    
    
    %
    c = Contrasts(1,:);
    cname = ContrastsNames{1};
    SPM.xCon = spm_FcUtil('Set', cname, 'T','c', c(:), SPM.xX.xKXs);
    
    for i=2:length(ContrastsNames)
        c = Contrasts(i,:);
        cname = ContrastsNames{i};
        SPM.xCon(i) = spm_FcUtil('Set', cname, 'T','c', c(:), SPM.xX.xKXs);
    end
    
    spm_contrasts(SPM);
    
    clear SPM
    
    cd (StartFolder)
    
    fprintf('\nThe analysis of the subject %s is done.\n\n', SubjID);
    
end