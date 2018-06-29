%%
StartFolder=fullfile(pwd, '..', '..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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


%%
for SubjInd = 1:size(SubjectList,1)
    
    clear ROI
    
    ROI(1) = struct('name', 'V1_surf_thres', 'fname', 'V1_surf_thres.nii','size', [], 'Tab', []);
    ROI(end+1) = struct('name', 'V2-3_surf_thres', 'fname', 'V2-3_surf_thres.nii','size', [], 'Tab', []);
    ROI(end+1) = struct('name', 'A1_surf', 'fname', 'A1_surf.nii','size', [], 'Tab', []);
    ROI(end+1) = struct('name', 'PT_surf_thres', 'fname', 'PT_surf_thres.nii','size', [], 'Tab', []);
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    Data_Folder = fullfile(SubjectFolder, 'Transfer');
    GlobalMask = spm_vol(fullfile(Data_Folder, 'rmask.nii'));
    GlobalMask = logical(spm_read_vols(GlobalMask));
    
    BackUpFolder = fullfile('/media','rxg243','BackUp2','AV_Integration_7T_2','Subjects_Data', ...
        ['Subject_' SubjID],'Transfer');
    Beta = dir(fullfile(BackUpFolder, 'r4beta*.nii'));
    Beta = spm_vol(fullfile(BackUpFolder, Beta(1).name));
    Beta = spm_read_vols(Beta);
    
    Struct_Folder = fullfile(SubjectFolder, 'Structural', 'CBS');
    Struct = dir(fullfile(Struct_Folder, 'T1*.nii'));
    Struct = spm_vol(fullfile(Struct_Folder, Struct.name));
    Struct = spm_read_vols(Struct);
    
    BrainCoverage(1) = sum(all([GlobalMask(:)==1 Struct(:)~=0],2))/sum(Struct(:)~=0); % whole brain coverage
    BrainCoverage(2) = sum(all([Beta(:)==0 Struct(:)~=0],2))/sum(Struct(:)~=0); % Number of empty values in beta file
    BrainCoverage(3) = sum(all([abs(Beta(:))<0.001 Beta(:)~=0 Struct(:)~=0],2))/sum(Struct(:)~=0); % Number of very small values in beta file
    BrainCoverage(4) = sum(all([abs(Beta(:))>0.001 Struct(:)~=0],2))/sum(Struct(:)~=0); % Number of normal values in beta file
    
    
    ROI_Folder = fullfile(Data_Folder, 'ROI');
    for iROI = 1:numel(ROI)
        RoiMask = spm_read_vols(spm_vol(fullfile(ROI_Folder,ROI(iROI).fname)));
        
        LeftROI = zeros(size(RoiMask));
        LeftROI(1:size(RoiMask,1)/2,:,:) = RoiMask(1:size(RoiMask,1)/2,:,:);
        
        ROI(iROI).size(1) = numel(find(LeftROI));
        
        ROI(iROI).Tab(1,1) = sum(Beta(find(LeftROI))==0)/sum(LeftROI(:));
        ROI(iROI).Tab(2,1) = (sum(abs(Beta(find(LeftROI)))<0.001) - sum(Beta(find(LeftROI))==0))/sum(LeftROI(:));
        ROI(iROI).Tab(3,1) = sum(abs(Beta(find(LeftROI)))>0.001)/sum(LeftROI(:));
        ROI(iROI).Tab(4,1) = sum(GlobalMask(find(LeftROI))==0)/sum(LeftROI(:));
        
        
        RightROI = zeros(size(RoiMask));
        RightROI(size(RoiMask,1)/2+1:end,:,:) = RoiMask(size(RoiMask,1)/2+1:end,:,:);
        
        ROI(iROI).size(2) = numel(find(RightROI));
        
        ROI(iROI).Tab(1,2) = sum(Beta(find(RightROI))==0)/sum(RightROI(:));
        ROI(iROI).Tab(2,2) = (sum(abs(Beta(find(RightROI)))<0.001) - sum(Beta(find(RightROI))==0))/sum(RightROI(:));
        ROI(iROI).Tab(3,2) = sum(abs(Beta(find(RightROI)))>0.001)/sum(RightROI(:));
        ROI(iROI).Tab(4,2) = sum(GlobalMask(find(RightROI))==0)/sum(RightROI(:));
        
    end
    
    save(fullfile(Data_Folder, 'ROI',['Subj_' SubjID '_ROI_coverage.mat']),'ROI')
    
    
    
end
