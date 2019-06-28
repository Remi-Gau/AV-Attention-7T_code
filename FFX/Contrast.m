% Computes some basic contrasts for some simple conditions and some basic differential contrasts 

clc; clear;

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

% StartFolder = fullfile(pwd, '..','..');
% StartFolder = '/media/rxg243/BackUp2/AV_Integration_7T_2';

StartFolder = '/media/remi/BackUp2/AV_Integration_7T_2';
CodeFolder = '/media/remi/BackUp2/AV_Integration_7T_2';
addpath(fullfile(CodeFolder, 'SubFun'))

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
    
    AnalysisFolder = fullfile(SubjectFolder, 'FFX');
    %     AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block');
    %     AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block_Smooth');
    %     AnalysisFolder = fullfile(SubjectFolder, 'FFX_NoSmooth');

    
    load(fullfile(AnalysisFolder, 'SPM.mat'))
    
    tmp = char(SPM.xX.name');
    tmp(:,1:6)=[];
    B=cellstr(tmp); clear tmp
    
    
    % SIMPLE CONTRASTS
    TEMP = zeros(1, size(SPM.xX.X,2));
    for CondInd=1:length(Conditions_Names)
        
        Contrasts = [Contrasts ; strcmp([Conditions_Names{CondInd} '*bf(1)'], B)']; %#ok<*AGROW>
        ContrastsNames{end+1}=strcat(Conditions_Names{CondInd}, ' > Baseline');  %#ok<*SAGROW>

    end

    TimeDer = ~cellfun('isempty',strfind(B,'*bf(2)'));
    
    AStim = ~cellfun('isempty',strfind(B,'A Stim'));
    AStim(TimeDer) = 0;
    
    AVStim = ~cellfun('isempty',strfind(B,'AV Stim'));
    AVStim(TimeDer) = 0;
    
    VStim = ~cellfun('isempty',strfind(B,'V Stim'));
    VStim(AVStim) = 0;
    VStim(TimeDer) = 0;
    
    
    ContrastsNames{end+1}='A_Stim > Baseline';  %#ok<*SAGROW>
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AStim) = 1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='AV Stim > Baseline';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AVStim) = 1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='V_Stim > Baseline';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(VStim) = 1;
    Contrasts = [Contrasts ; TEMP];
    


    % DIFFERENTIAL CONTRASTS
    ContrastsNames{end+1}='A Stim - A Att > A Stim - V Att';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(strcmp('A Stim - Auditory Attention*bf(1)', B)) = 1;
    TEMP(strcmp('A Stim - Visual Attention*bf(1)', B)) = -1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='V Stim - A Att > V Stim - V Att';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(strcmp('V Stim - Auditory Attention*bf(1)', B)) = 1;
    TEMP(strcmp('V Stim - Visual Attention*bf(1)', B)) = -1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='AV Stim - A Att > AV Stim - V Att';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(strcmp('AV Stim - Auditory Attention*bf(1)', B)) = 1;
    TEMP(strcmp('AV Stim - Visual Attention*bf(1)', B)) = -1;
    Contrasts = [Contrasts ; TEMP];
    

    ContrastsNames{end+1}='AV > A + V';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AVStim) = 1;
    TEMP(VStim) = -1;
    TEMP(AStim) = -1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='AV - A > 0';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AVStim) = 1;
    TEMP(AStim) = -1;
    Contrasts = [Contrasts ; TEMP];
    
    ContrastsNames{end+1}='AV - V > 0';
    TEMP = zeros(1, size(SPM.xX.X,2));
    TEMP(AVStim) = 1;
    TEMP(VStim) = -1;
    Contrasts = [Contrasts ; TEMP];
    

    % set the batch
    for i=1:length(ContrastsNames)
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(AnalysisFolder, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = ContrastsNames{i};
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = double(Contrasts(i,:));
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 1;
    end
    
    spm_jobman('run', matlabbatch)
    
    fprintf('\nThe analysis of the subject %s is done.\n\n', SubjID);
    
end