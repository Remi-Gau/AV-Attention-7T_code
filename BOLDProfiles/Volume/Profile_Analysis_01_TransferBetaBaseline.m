clear; clc;

% StartDirectory = fullfile(pwd, '..', '..');
StartDirectory = '/media/rxg243/BackUp2/AV_Integration_7T_2';

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

CondNames = {...
    'A Stim - Auditory Attention';....
    'V Stim - Auditory Attention';....
    'AV Stim - Auditory Attention';....
    'A Stim - Visual Attention';....
    'V Stim - Visual Attention';....
    'AV Stim - Visual Attention';....
    'Baseline'}


for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    mkdir(fullfile(SubjectFolder, 'FFX_baseline', 'Beta'))
    
    load(fullfile(SubjectFolder, 'FFX_baseline', 'SPM.mat'))
    
    BetaNames = char(SPM.xX.name');
    
    BetaOfInterest = [];
    
    for iCond = 1:numel(CondNames)
        tmp = strfind(cellstr(BetaNames),[CondNames{iCond} '*bf(1)']);
        BetaOfInterest(:,iCond) = ~cellfun('isempty',tmp);
        clear tmp
    end
    
    BetaOfInterest = find(any(BetaOfInterest,2));
    
    %     for i=1:length(BetaOfInterest)
    %         copyfile(fullfile(SubjectFolder, 'FFX_baseline',  sprintf('beta_%04d.nii', BetaOfInterest(i))), ...
    %             fullfile(SubjectFolder, 'FFX_baseline', 'Beta'))
    %     end
    %
    %     copyfile(fullfile(SubjectFolder, 'FFX_baseline', 'mask.nii' ), ...
    %             fullfile(SubjectFolder, 'FFX_baseline', 'Beta'))
    
    %         copyfile(fullfile(SubjectFolder, 'mean*.nii'), fullfile(SubjectFolder, 'FFX_baseline', 'Beta'))
    
    copyfile(fullfile(SubjectFolder, 'Nifti', 'NoMoCo', '01', 'meanUR*.nii'), fullfile(SubjectFolder, 'FFX_baseline', 'Beta'))
end