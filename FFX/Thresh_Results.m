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
    'A Stim', ...
    'A Stim', ...
    'V_AV Stim',...
    'V_AV Stim',...
    'V Stim', ...
    'V Stim'};


for SubjInd = 1:size(SubjectList,1)
    
    ContrastsNames={};
    Contrasts = [];
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    AnalysisFolder = fullfile(SubjectFolder, 'FFX_Block_Smooth');
    
    for iCond = 1:6
        if mod(iCond,2)~=0
            tmp = '>Baseline';
        else
            tmp = '<Baseline';
        end
        
        for Thresh = 0.1
            
            matlabbatch{1}.spm.stats.results.spmmat = {fullfile(AnalysisFolder,'SPM.mat')};
            matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
            matlabbatch{1}.spm.stats.results.conspec.contrasts = iCond;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.print = false;
            matlabbatch{1}.spm.stats.results.write.binary.basename = ...
                [Conditions_Names{iCond} tmp '_thres_' num2str(Thresh)];
            
            spm_jobman('run', matlabbatch)
            
            fprintf('\nThe analysis of the subject %s is done.\n\n', SubjID);
            
        end
    end
    
end
