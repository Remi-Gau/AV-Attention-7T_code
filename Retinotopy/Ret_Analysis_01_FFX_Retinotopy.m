%%
clear all; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>
defaults.mask.thresh    = -1;

PrefixImages2Select = 'UR*.nii';

CondNamesTotal = {'Sine', 'Cosine'};

% Time Slicing Parameters
TR = 3;
NbSlices = 60;
HPF = 128;
ReferenceSlice = 1;

% Retinotopy parameters
%NbVolPerCycle = 18.3; % Check PTB scripts for that (usually NbVolPerCycle * TR ~ 1 min)
NbCycles=8;

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

%  Root folder definition
StartFolder = pwd;


%% Define folders, number of runs, scans per run...
for SubjInd = 1:size(SubjectList,1)
    
    NbVol2Discard = [];
    NbVol = [];
    NbVolPerCycle = [];
    Direction = [];
    
    tic
    
    %% Subject's Identity and folders
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiFolder = fullfile(SubjectFolder, 'Retinotopy', 'Nifti');
    
    AnalysisFolder = fullfile(SubjectFolder, 'Retinotopy', 'FFX');
    
    copyfile(fullfile(AnalysisFolder, 'mask.nii'), fullfile(SubjectFolder, 'Retinotopy'))
    
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
        
        Stim_Time={};
        
        Stim_Time{1,1} = FileContent{1,3};
        Stim_Time{1,2} = char(FileContent{1,4});
        
        % Number of volumes to discard and total number of volumes to keep
        Triggers = find(strcmp(Stim_Time{1,1}, '30'));
        
        Start = find(strcmp(Stim_Time{1,1}, 'NewCycle'), 1, 'first');
        
        NbVol2Discard(iFile) = sum(Triggers<Start)-1;
        StartXP = str2num(Stim_Time{1,2}(Triggers(NbVol2Discard(iFile)),:));
        StartTime(iFile) = (str2num(Stim_Time{1,2}(Start,:))-StartXP)/10000;
        
        NbVol(iFile) = numel(Triggers)-NbVol2Discard(iFile);
        
        DurationRun(iFile) = (str2num(Stim_Time{1,2}(strcmp('Final_Fixation', Stim_Time{1,1}),:))-StartXP)/10000;
        
        % Shortens realignment parameters file accordingly
        cd(fullfile(NiftiFolder, sprintf('%2.2d', iFile)))
        Mov_Parameter_ls = dir('rp_*.txt');
        R = load(fullfile(pwd, Mov_Parameter_ls.name));
        R = R(NbVol2Discard(iFile)+(1:NbVol(iFile)),:);
        save('rp.mat', 'R')
        clear R Mov_Parameter_ls
        
        % Number of volumes per cycle
        TEMP = find(strcmp(Stim_Time{1,1}, '30'));
        TEMP2 = find(strcmp(Stim_Time{1,1}, 'NewCycle'));
        
        for i=1:length(TEMP2)-1
            TEMP3(i) = sum(all([TEMP>TEMP2(i) TEMP<TEMP2(i+1)],2));
        end
        NbVolPerCycle(iFile) = mean(TEMP3);
        clear TEMP3
        
        % Analyse accuracy and RT
        %         Stim_Time{1,1}(1:Start) = [];
        %         Stim_Time{1,1}(1:Start) = [];
        %
        %         Stim_Time{1,1}([TEMP ; TEMP2]) = [];
        %         Stim_Time{1,2}([TEMP ; TEMP2]) = [];
        
    end
    
    %% Specify the batch
    mkdir(AnalysisFolder)
    cd(AnalysisFolder)
    delete *.*
    
    fprintf('\nSpecifying the job\n\n')
    
    matlabbatch ={};
    
    matlabbatch{1,1}.spm.stats.fmri_spec.dir{1,1} = AnalysisFolder;
    
    matlabbatch{1,1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1,1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1,1}.spm.stats.fmri_spec.timing.fmri_t = NbSlices;
    matlabbatch{1,1}.spm.stats.fmri_spec.timing.fmri_t0 = ReferenceSlice;
    
    matlabbatch{1,1}.spm.stats.fmri_spec.fact = struct('name',{},'levels',{});
    
    matlabbatch{1,1}.spm.stats.fmri_spec.bases.hrf.derivs = [0,0]; % First is time derivative, Second is dispersion
    
    matlabbatch{1,1}.spm.stats.fmri_spec.volt = 1;
    
    matlabbatch{1,1}.spm.stats.fmri_spec.global = 'None';
    
    matlabbatch{1,1}.spm.stats.fmri_spec.mask = {fullfile(SubjectFolder, 'Retinotopy', 'mask.nii')};
    
    matlabbatch{1,1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    for RunInd = 1:length(NbVol)
        
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).multi{1} = '';
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).regress = struct('name',{},'val',{});
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).hpf = HPF;
        
        % Enter source folder reads the image files
        cd(fullfile(NiftiFolder, sprintf('%2.2d', RunInd)))
        
        % Lists the images
        IMAGES_ls = dir(PrefixImages2Select);
        IMAGES_ls = spm_vol(IMAGES_ls.name);
        % Names them with their absolute pathnames
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).scans = {};
        for j = NbVol2Discard(RunInd)+(1:NbVol(RunInd))
            matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).scans{end+1,1} = [fullfile(pwd, filesep, IMAGES_ls(j).fname), ',' , num2str(j)];
        end
        
        % Names them with its absolute pathname
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).multi_reg{1} = fullfile(pwd, 'rp.mat');
        
        % Onsets
        onset_vec = StartTime(RunInd):TR:DurationRun(RunInd);
        regressor{1} = sin(2*pi*1/(NbVolPerCycle(RunInd)*TR).*onset_vec);
        regressor{2} = cos(2*pi*1/(NbVolPerCycle(RunInd)*TR).*onset_vec);
        
        % Conditions
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.name = 'Onsets';
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.onset = onset_vec;
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.duration = 0;
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.tmod = 0;
        
        % Parametric modulation
        for ParamModInd = 1:length(CondNamesTotal)
            
            if Direction(RunInd)=='+'
                order = 2;
            else
                order = 1;
            end
            
            matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.pmod(ParamModInd).name = CondNamesTotal{ParamModInd};
            
            if ParamModInd == 1
                if order == 2
                    matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.pmod(ParamModInd).param = -regressor{ParamModInd}';
                else
                    matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.pmod(ParamModInd).param = regressor{ParamModInd}';
                end
            else
                matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.pmod(ParamModInd).param = regressor{ParamModInd}';
            end
            
            matlabbatch{1,1}.spm.stats.fmri_spec.sess(1,RunInd).cond.pmod(ParamModInd).poly = 1;
        end
        
        clear onset_vec regressor
    end
    
    
    %% fMRI estimation
    fprintf('\nSpecifying and estimating model\n\n')
    matlabbatch{1,end+1}={};
    matlabbatch{1,end}.spm.stats.fmri_est.spmmat{1,1} = fullfile(AnalysisFolder, 'SPM.mat');     %set the spm file to be estimated
    matlabbatch{1,end}.spm.stats.fmri_est.method.Classical = 1;
    
    cd(AnalysisFolder)
    save (strcat('FFX_Retinotopy_Subject_', SubjID, '_jobs'));
    
    spm_jobman('run', matlabbatch)
    fprintf('\nThe analysis of the subject %s is done.\n\n', SubjID);
    
    cd (StartFolder)
    
    toc
    
end
