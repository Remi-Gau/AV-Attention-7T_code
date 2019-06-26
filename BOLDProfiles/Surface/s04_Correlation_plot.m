% plots the vertex wise correlation 
% fisher transformation applied
% significance computed with exact sign permutation test

%%
clc; clear; close all

CodeFolder = '/home/remi/github/AV-Attention-7T_code';
addpath(genpath(fullfile(CodeFolder, 'SubFun')))
Get_dependencies('/home/remi')

SourceFolder = '/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T';

opt.FigDim = [100, 100, 1000, 1500];
opt.Visibility = 'on';
opt.Xpos =  [1 1.5];
opt.Subj2Include = true(11,1);
opt.fontsize = 10;
opt.FigureFolder = fullfile(CodeFolder,'Figures');
opt.print = 0;

load(fullfile(SourceFolder,'Results','Profiles','Surfaces','SurfCorrelation.mat'))

%% Basic conditions
CdtionName={'A vs V', 'A vs AV', 'V vs AV'};
opt.name = 'Grp_CorrBaseline';

plot_group_vertex_correlation(GrpBetaCdt, ROI, CdtionName, opt)


%% Cross-sens
CdtionName={'A vs AV-V','V vs AV-A'};
opt.name = 'Grp_CrossSensory';

plot_group_vertex_correlation(GrpBetaCrossSens, ROI, CdtionName, opt)


%% Attention
CdtionName={...
    'Att Mod_{A stim} vs Att Mod_{V stim}', ...
    'Att Mod_{A stim} vs Att Mod_{AV stim}', ...
    'Att Mod_{V stim} vs AttMod_{AV stim}'};
opt.name = 'Grp_CorrAttention';

plot_group_vertex_correlation(GrpBetaAtt, ROI, CdtionName, opt)


%% Stim vs Att
CdtionName={'A stim vs Att Mod','V stim vs Att Mod','AV Stim vs AttMod'};
opt.name = 'Grp_Corr_Stim_VS_Attention';

plot_group_vertex_correlation(GrpBetaStimAtt, ROI, CdtionName, opt)


%% Stim vs stim specific Att
CdtionName={'A stim vs Att Mod A stim','V stim vs Att Mod V stim','AV Stim vs AttMod AV stim'};
opt.name = 'Grp_Corr_StimSpe_VS_Attention';

plot_group_vertex_correlation(GrpBetaStimSpeAtt, ROI, CdtionName, opt)


%% Cross vs Att
CdtionName={'AV-A vs Att Mod','AV-V vs Att Mod'};
opt.name = 'Grp_Corr_CrossMod_VS_Attention';

plot_group_vertex_correlation(GrpCrossModAtt, ROI, CdtionName, opt)
