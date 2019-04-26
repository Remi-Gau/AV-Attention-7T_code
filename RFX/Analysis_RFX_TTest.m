% Do group level analysis (one sample ttest) of the basic contrasts
% computed at the subject level

StartFolder = '/media/remi/BackUp2/AV_Integration_7T_2';
CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(fullfile(CodeFolder, 'SubFun'))

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


%% create mean structural

for SubjInd = 1:size(SubjectList,1)
    Struct{SubjInd,1} = fullfile(StartFolder, 'Subjects_Data', ...
        ['Subject_' SubjectList(SubjInd,:)], ...
        'Structural', 'SPM_Segmentation', 'wmUNI.nii');
end


temp = spm_vol(char(Struct));
vol = mean(spm_read_vols(temp),4);
hdr = temp(1);
hdr.fname = fullfile(StartFolder, 'RFX', 'AvgStruct.nii');
spm_write_vol(hdr, vol);


%%
matlabbatch = {};

load(fullfile(StartFolder, 'Subjects_Data/Subject_02/FFX/SPM.mat'))

for iTest = 1:numel(SPM.xCon)
    
    ContrastsNames = SPM.xCon(iTest).name
    
    mkdir(fullfile(StartFolder, 'RFX', ContrastsNames));
    delete(fullfile(StartFolder, 'RFX', ContrastsNames, 'SPM.mat'));

    % Design
    matlabbatch{end+1}.spm.stats.factorial_design.dir = ...
        {fullfile(StartFolder, 'RFX', ContrastsNames)};
    
    matlabbatch{end}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{end}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{end}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{end}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{end}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{end}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{end}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{end}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    for SubjInd = 1:size(SubjectList,1)
        matlabbatch{end}.spm.stats.factorial_design.des.t1.scans{SubjInd,1} = ...
            fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjectList(SubjInd,:)], ...
            'FFX', sprintf('con_%4.4i.nii', iTest));
    end
    
    % Estimate
    matlabbatch{end+1}={};
    matlabbatch{end}.spm.stats.fmri_est.spmmat{1} = fullfile(StartFolder, ...
        'RFX', ContrastsNames, 'SPM.mat');
    matlabbatch{end}.spm.stats.fmri_est.method.Classical = 1;
    
    % Contrast
    matlabbatch{end+1}.spm.stats.con.spmmat{1} = fullfile(StartFolder, ...
        'RFX', ContrastsNames, 'SPM.mat');
    matlabbatch{end}.spm.stats.con.consess{1}.tcon.name = ContrastsNames;
    matlabbatch{end}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{end}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{end}.spm.stats.con.delete = 1;
    
    % NIDM results
    matlabbatch{end+1}.spm.stats.results.spmmat{1} = fullfile(StartFolder, ...
        'RFX', ContrastsNames, 'SPM.mat');
    matlabbatch{end}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{end}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{end}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{end}.spm.stats.results.conspec.thresh = 0.005;
    matlabbatch{end}.spm.stats.results.conspec.extent = 10;
    matlabbatch{end}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{end}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{end}.spm.stats.results.units = 1;
    matlabbatch{end}.spm.stats.results.export{1}.ps = true;
    matlabbatch{end}.spm.stats.results.export{2}.nidm.modality = 'FMRI';
    matlabbatch{end}.spm.stats.results.export{2}.nidm.refspace = 'ixi';
    matlabbatch{end}.spm.stats.results.export{2}.nidm.group.nsubj = 11;
    matlabbatch{end}.spm.stats.results.export{2}.nidm.group.label = 'ctrl';
    
end

spm_jobman('run', matlabbatch)
