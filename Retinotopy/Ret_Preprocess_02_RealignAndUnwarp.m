% Add an option for the Fielmap correction
clear all; clc

spm_jobman('initcfg')
spm_get_defaults;
global defaults %#ok<NUSED>


%  Folders definitions
RootFolder = pwd;

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
    '16';...
    ];


FoldersNames = {...
    1;...
    1:2;...
    1:2;...
    1;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    1:2;...
    };

% if ~matlabpool('SIZE')
%     matlabpool('OPEN', length(SubjectList)) ;
% end


for SubjInd = 7%:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:)
    
    SubjectFolder = fullfile(pwd, 'Subjects_Data', ['Subject_' SubjID]);
    
    NiftiSourceFolder = fullfile(SubjectFolder, 'Retinotopy', 'Nifti');
    
    FieldmapsFolder = fullfile(SubjectFolder, 'Retinotopy', 'FieldMap');
    
    NbRuns = length( FoldersNames{SubjInd});
    
    
    %% ---------------------------  %
    %      UNWARP & REALIGN         %
    %  ---------------------------  %
    
    tic
    
    matlabbatch = {};
    
    for RunInd=1:NbRuns
        
        RunInd
        
        cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(RunInd))))
        TEMP = dir('S*iPAT3_6_8_60sli_TE18ms_1p5_wb_Te18.nii');
        TEMP = spm_vol(TEMP.name);
        
        for i=1:length(TEMP)
            matlabbatch{1,1}.spm.spatial.realignunwarp.data(1,RunInd).scans{i} = [fullfile(pwd,TEMP(i).fname) ,',', num2str(i)];
        end
        
        if SubjID == '11'
            matlabbatch{1,1}.spm.spatial.realignunwarp.data(1,RunInd).pmscan = {''};
        else
            cd(FieldmapsFolder)
            TEMP=dir('vdm*.nii')
            matlabbatch{1,1}.spm.spatial.realignunwarp.data(1,RunInd).pmscan = {[fullfile(pwd,TEMP.name) ',1']};
        end
    end
    
    % --------------------------%
    %     DEFINES    BATCH      %
    % --------------------------%
    
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.quality = 1;
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.sep = 1.5; % Defaut = 4
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.fwhm = 1; % Default = 5
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.rtm = 0; % Register to mean: 0 ---> FALSE
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 1 0];
    matlabbatch{1,1}.spm.spatial.realignunwarp.eoptions.weight = {''};
    
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    
    matlabbatch{1,1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1,1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 1 0];
    matlabbatch{1,1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1,1}.spm.spatial.realignunwarp.uwroptions.prefix = 'UR';
    
    cd(fullfile(NiftiSourceFolder, sprintf('%2.2d', FoldersNames{SubjInd}(1))))
    save (strcat('RealAndUnwarpWithFieldMap_Subj_', SubjID, '_matlabbatch'));
    
    fprintf('\n\n')
    disp('%%%%%%%%%%%%%%%%')
    disp('   REALIGNING   ')
    disp('%%%%%%%%%%%%%%%%')
    
    spm_jobman('run',matlabbatch)
    
    cd (RootFolder)
    
    toc
    
end

% matlabpool('CLOSE')