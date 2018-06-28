%%
clear all; clc

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
   
    AnalysisFolder = fullfile(SubjectFolder, 'Retinotopy', 'FFX');

    %%
    load(fullfile(AnalysisFolder, 'SPM.mat'))
    
    no_of_runs = length(SPM.Sess);   
    Polar_sine_columns = [];
    Polar_cosine_columns = [];
    Polar_BETA_images = cell(no_of_runs,2);
    
    for iRun = 1:no_of_runs
            columns_sess =  SPM.Sess(iRun).col; %#ok<SAGROW>
            disp(SPM.xX.name{1,columns_sess(2)})
            Polar_sine_columns = [Polar_sine_columns columns_sess(2)]; %#ok<AGROW>
            disp(SPM.xX.name{1,columns_sess(3)})
            Polar_cosine_columns = [Polar_cosine_columns columns_sess(3)]; %#ok<AGROW>
    end
    
    for iRun = 1:no_of_runs
        Polar_BETA_images{iRun,1} = fullfile(AnalysisFolder,sprintf('beta_%04d.nii', Polar_sine_columns(iRun)));
        Polar_BETA_images{iRun,2} = fullfile(AnalysisFolder,sprintf('beta_%04d.nii', Polar_cosine_columns(iRun)));
    end
    
    Mask = logical(spm_read_vols(spm_vol(fullfile(AnalysisFolder, 'mask.nii'))));
    
    %% POLAR
    % reading headers for images of interest
    Sine_hdr = spm_vol(str2mat(Polar_BETA_images{:,1})); %#ok<FPARK>
    Cosine_hdr = spm_vol(str2mat(Polar_BETA_images{:,2})); %#ok<FPARK>

    % read volumes of the respective headers
    Sine_vols = spm_read_vols(Sine_hdr);
    Cosine_vols = spm_read_vols(Cosine_hdr);
    
    % average beta images for of sine/cosine regressors 
    % (note: clock/anticlockwise directions have been accounted for in the
    % FFX model so no need to care about them here)
    mean_SineBetas = mean(Sine_vols,4);
    mean_CosineBetas = mean(Cosine_vols,4);
    
    clear Sine_vols Cosine_vols
    
    % compute phase and amplitude
    Y_phase = mod(atan2(mean_SineBetas, mean_CosineBetas)*180/pi, 360);
    Y_phase(~Mask) = NaN;
    Y_ampl = sqrt(mean_SineBetas.^2 + mean_CosineBetas.^2);
    Y_ampl(~Mask) = NaN;
       
    mean_SineBetas(~Mask) = NaN;
    mean_SineBetas(~Mask) = NaN;
    
    % taking header information from one of the images above
    newImgInfo = Sine_hdr(1);

    % writing image for mean beta for sine regressor
    meanSineBetas_hdr = newImgInfo;
    meanSineBetas_hdr.fname = fullfile(AnalysisFolder, 'Polar_MeanSineBetas.nii');
    meanSineBetas_hdr.private.dat.fname = meanSineBetas_hdr.fname;
    spm_write_vol(meanSineBetas_hdr,mean_SineBetas);

    % writing image for mean beta for sine regressor
    meanCosineBetas_hdr = newImgInfo;
    meanCosineBetas_hdr.fname = fullfile(AnalysisFolder, 'Polar_MeanCosineBetas.nii');
    meanCosineBetas_hdr.private.dat.fname = meanCosineBetas_hdr.fname;
    spm_write_vol(meanCosineBetas_hdr,mean_CosineBetas);
    
    % writing image for phase
    phase_hdr = newImgInfo;
    phase_hdr.fname = fullfile(AnalysisFolder, 'Polar_Phase.nii');
    phase_hdr.private.dat.fname = phase_hdr.fname;
    spm_write_vol(phase_hdr, Y_phase);
    
    % writing image for amplitude
    amplitude_hdr = newImgInfo;
    amplitude_hdr.fname = fullfile(AnalysisFolder, 'Polar_Amplitude.nii');
    amplitude_hdr.private.dat.fname = amplitude_hdr.fname;
    spm_write_vol(amplitude_hdr, Y_ampl);
    
    clear Y_ampl Y_phase mean_SineBetas mean_CosineBetas 

    
end


