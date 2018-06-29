clear; clc

StartFolder = pwd;

SubjectList = [...
    '02';...
%     '03';...
    '04';...
    %'06';...
    %'07';...
%     '08';...
%     '09';...
%     '11';...
    '12';...
%     '13';...
    %'14';...
%     '15';...
%     '16'
    ];


NbWorkers = 3;
if isempty(gcp)
    parpool(NbWorkers)
end

parfor SubjInd = 1:size(SubjectList,1)
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    fprintf('\n\nAnalysing subject %s\n', num2str(SubjID))
    
    %%  Structural
    Structural = fullfile(SubjectFolder, 'Structural', 'CBS', 'Layering', 'T1_03_Layers.nii');
    LayerLabels = spm_read_vols(spm_vol(Structural));
    Layers = unique(LayerLabels(:));
    Layers(Layers==0) = [];
    
        
    %% Beta 
    
    AnalysisFolder = fullfile(SubjectFolder, 'Transfer');
    
    cd(AnalysisFolder)
    
    FileList = dir('rbeta*.nii');
    
    for iFile = 1:numel(FileList)
        fprintf('  %s\n', FileList(iFile).name)
        
        Hdr = spm_vol(FileList(iFile).name);
        Vol = spm_read_vols(Hdr);
        
        for FWHM = [3 6];
            
            VolFinal= nan(size(Vol));
            
            for iLayer= 1:numel(Layers)
                
                VolTmp = nan(size(Vol));
                VolTmp(ismember(LayerLabels,Layers(iLayer))) = ...
                    Vol(ismember(LayerLabels,Layers(iLayer)));
                
                spm_smooth(VolTmp,VolTmp,FWHM,0)

                VolFinal(ismember(LayerLabels,Layers(iLayer))) = ...
                    VolTmp(ismember(LayerLabels,Layers(iLayer)));

            end
            
            HdrTmp = Hdr;
            HdrTmp.fname = ['S' num2str(FWHM) HdrTmp.fname];
            
            spm_write_vol(HdrTmp, VolFinal);
            
        end
        
    end
    
    
    cd(StartFolder)
end