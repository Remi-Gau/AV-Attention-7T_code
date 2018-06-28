%%
clear all; clc; close all;

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

%  Root folder definition
StartFolder = pwd;

cd(StartFolder)

for SubjInd = 1:size(SubjectList,1)
    %% Subject's Identity and folders
    SubjID = SubjectList(SubjInd,:) %#ok<NOPTS>
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
   
    AnalysisFolder = fullfile(SubjectFolder, 'Retinotopy', 'MIPAV');

   
    %% POLAR
    % reading headers for images of interest
    Cosine_hdr = spm_vol(fullfile(AnalysisFolder, 'rPolar_MeanCosineBetas.nii')); %#ok<FPARK>
    Sine_hdr = spm_vol(fullfile(AnalysisFolder, 'rPolar_MeanSineBetas.nii')); %#ok<FPARK>

    % read volumes of the respective headers
    Sine_vols = spm_read_vols(Sine_hdr);
    Cosine_vols = spm_read_vols(Cosine_hdr);
    
    % compute phase and amplitude
    Y_phase = atan2(Sine_vols, Cosine_vols);
    Y_ampl = sqrt(Sine_vols.^2 + Cosine_vols.^2);
    
    % taking header information from one of the images above
    newImgInfo = Sine_hdr;
    
    % writing image for phase
    phase_hdr = newImgInfo;
    phase_hdr.fname = fullfile(AnalysisFolder, 'rPolar_Phase.nii');
    phase_hdr.private.dat.fname = phase_hdr.fname;
    spm_write_vol(phase_hdr, Y_phase);
    
    % writing image for amplitude
    amplitude_hdr = newImgInfo;
    amplitude_hdr.fname = fullfile(AnalysisFolder, 'rPolar_Amplitude.nii');
    amplitude_hdr.private.dat.fname = amplitude_hdr.fname;
    spm_write_vol(amplitude_hdr, Y_ampl);
    
    clear Y_ampl Y_phase mean_SineBetas mean_CosineBetas 

    
end


