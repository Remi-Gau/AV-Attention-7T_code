clear;
clc;

StartDirectory = pwd;

%%
NbLayer=6;
Layer2Display = {6; 5; 4; 3; 2; 1};
% Layer2Display = {[6 5];[4];[3 2 1]};
% Layer2Display = {[6 5];[4 3];[2 1]};

% NbLayer=10;
% Layer2Display = {10; 9; 8; 7; 6; 5; 4; 3; 2; 1};

SubjID = 1;
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

SubROI_Nb = 19;
SubROI_List = {...
    'V1L' , 'V1_L_MNI.nii', 'V_lh'; ... % 1
    'V2L' , 'V2_L_MNI.nii', 'V_lh'; ...
    'V3dL' , 'V3d_L_MNI.nii', 'V_lh'; ...
    'V3vL' , 'V3v_L_MNI.nii', 'V_lh'; ...
    'V4dL' , 'V4d_L_MNI.nii', 'V_lh'; ... % 5
    'V4vL' , 'V4v_L_MNI.nii', 'V_lh'; ...
    'V5L' , 'V5_L_MNI.nii', 'V_lh'; ...
    
    'V1R' , 'V1_R_MNI.nii', 'V_rh'; ...
    'V2R' , 'V2_R_MNI.nii', 'V_rh'; ...
    'V3dR' , 'V3d_R_MNI.nii', 'V_rh'; ... % 10
    'V3vR' , 'V3v_R_MNI.nii', 'V_rh'; ...
    'V4dR' , 'V4d_R_MNI.nii', 'V_rh'; ...
    'V4vR' , 'V4v_R_MNI.nii', 'V_rh'; ...
    'V5vR' , 'V5_R_MNI.nii', 'V_rh'; ...
    
    'PTL' , 'PT_L_MNI.nii', 'A_lh'; ... % 15
    'TE1.0L' , 'TE_1.0_L_MNI.nii', 'A_lh'; ...
    'TE1.1L' , 'TE_1.1_L_MNI.nii', 'A_lh'; ...
    'TE1.2L' , 'TE_1.2_L_MNI.nii', 'A_lh'; ...
    'A1L' , 'A1_L.nii', 'A_lh'; ...
    
    'PTR' , 'PT_R_MNI.nii', 'A_rh'; ... % 20
    'TE1.0R' , 'TE_1.0_R_MNI.nii', 'A_rh'; ...
    'TE1.1R' , 'TE_1.1_R_MNI.nii', 'A_rh'; ...
    'TE1.2R' , 'TE_1.2_R_MNI.nii', 'A_rh'; ...
    'A1R' , 'A1_R.nii', 'A_rh'};


Analysis = struct('name', 'A Stim VS V Stim');
Analysis(end+1) = struct('name', 'A Stim VS AV Stim');
Analysis(end+1) = struct('name', 'V Stim VS AV Stim');
Analysis(end+1) = struct('name', 'A Att VS V Att');


AnalysisFolder=fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjectList(SubjID,:)], ...
    'Analysis', 'ROI', SubROI_List{SubROI_Nb,3});


%%
%[BlobImg, sts] = spm_select([1 1],'image','Select contrast');
BlobImg = fullfile(AnalysisFolder, 'SVM', [SubROI_List{SubROI_Nb,1} '_weight_1.nii']);
if ischar(BlobImg), BlobHdr = spm_vol(BlobImg);  BlobVol = spm_read_vols(BlobHdr); end

%[MaskImg, sts] = spm_select([1 1],'image','Select mask');
MaskImg = fullfile(AnalysisFolder, ['rsub_' SubROI_List{SubROI_Nb,2}]);
if ischar(MaskImg), MaskHdr = spm_vol(MaskImg);  MaskVol = spm_read_vols(MaskHdr); end

BlobVol(~MaskVol)=0;
BlobVol(isnan(BlobVol))=0;

%[MaskImg, sts] = spm_select([1 1],'image','Select mask');
LayerImg = fullfile(AnalysisFolder, ['T1_' sprintf('%0.2d',NbLayer) '_Layers.nii']);
if ischar(LayerImg), LayerHdr = spm_vol(LayerImg);  LayerVol = spm_read_vols(LayerHdr); end

LayerVol(~MaskVol)=0;

%[images, sts] = spm_select([1 1],'image','Select structural');
images = fullfile(AnalysisFolder, 'T1.nii');
images = repmat(images, NbLayer, 1);
if ischar(images), images = spm_vol(images); end

%%
Fig = spm_figure('GetWin','Graphics');
WinSca = spm('WinScale');
spm_figure('Clear','Graphics');
spm_orthviews('Reset');

uicontrol(Fig,'style','text','string',['Layers overlay; Subject: ' SubjectList(SubjID,:) ...
    '; ROI: ' SubROI_List{SubROI_Nb,1}],'position', ...
    [0 850 600 15].*WinSca,'Fontsize',13,'backgroundcolor',[1 1 1]);

%%
n  = round(length(Layer2Display)^0.4);
m  = ceil(length(Layer2Display)/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;

for LayerGroupInd=1:length(Layer2Display)
    
    i = 1-h*(floor((LayerGroupInd-1)/n)+1);
    j = w*rem(LayerGroupInd-1,n);
    
    handle = spm_orthviews('Image', images(LayerGroupInd),...
        [j+ds/2 i+ds/2 w-ds h-ds]);
    
    I=[];
    for LayInd=1:length(Layer2Display{LayerGroupInd})
        I = [I ; find(LayerVol==Layer2Display{LayerGroupInd}(LayInd)) ];
    end;
    
    t=BlobVol(I);
    [X,Y,Z] = ind2sub(size(BlobVol),I);
    
    spm_orthviews('AddBlobs',handle, [X';Y';Z'], t', BlobHdr.mat, '');
    %spm_orthviews('Caption', handle, '                B values');
    
    j2=w*rem(LayerGroupInd,n);
    uicontrol(Fig,'style','text','string', ['Layers ' num2str(Layer2Display{LayerGroupInd})],'position', ...
        [600*j2 825-840*i 600/n 15].*WinSca,'Fontsize',13,'backgroundcolor',[1 1 1]);
    
    clear X Y Z t
    
end