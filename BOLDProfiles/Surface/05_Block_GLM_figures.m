%%
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

Smooth=1;

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

hs = 'lr';
Sens = {'AStim', 'VStim', 'AVStim'};
CrossSens = {'AV-A', 'AV-V'};
Att  = {'A_Stim_A_Att-V_Att', 'V_Stim_A_Att-V_Att', 'AV_Stim_A_Att-V_Att', 'A_Att-V_Att'};
% Att  = {'A_Stim_A_Att-V_Att', 'V_Stim_A_Att-V_Att', 'AV_Stim_A_Att-V_Att'};

Var={'Cst','Lin'};

if Smooth
    Sufix = '_smooth_6fwhm';
    CLIM = {[-1.5 1.5], [-.5 .5]};
else
    Sufix = ''; %#ok<*UNRCH>
    CLIM = {[-3 3], [-1 1]};
end

% Color map
X = 0:0.001:1;
R = 0.237 - 2.13*X + 26.92*X.^2 - 65.5*X.^3 + 63.5*X.^4 - 22.36*X.^5;
G = ((0.572 + 1.524*X - 1.811*X.^2)./(1 - 0.291*X + 0.1574*X.^2)).^2;
B = 1./(1.579 - 4.03*X + 12.92*X.^2 - 31.4*X.^3 + 48.6*X.^4 - 23.36*X.^5);

figure('name','Legend')
subplot(121)
colormap([R' G' B']);
imagesc(repmat([2:-.01:-2]', [1,200]), [-2 2])
title('Constant')

subplot(122)
colormap([R' G' B']);
imagesc(repmat([2:-.01:-2]', [1,200]), [-2 2])
title('Linear')

mkdir(fullfile(StartFolder,'Figures','Profiles', 'Surfaces'))

subplot(121)
set(gca,'xtick',[],'ytick',linspace(1,400,9),'yticklabel',linspace(3,-3,9));
subplot(122)
set(gca,'xtick',[],'ytick',linspace(1,400,9),'yticklabel',linspace(1,-1,9));
print(gcf, fullfile(StartFolder,'Figures','Profiles', 'Surfaces', 'Scales_CstLin_Surf.tif'), '-dtiff')

subplot(121)
set(gca,'xtick',[],'ytick',linspace(1,400,9),'yticklabel',linspace(1.5,-1.5,9));
subplot(122)
set(gca,'xtick',[],'ytick',linspace(1,400,9),'yticklabel',linspace(.5,-.5,9));
print(gcf, fullfile(StartFolder,'Figures','Profiles', 'Surfaces', 'Scales_CstLin_Surf_Smooth.tif'), '-dtiff')

%%
for SubjInd = 2:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    Results_Folder = fullfile(SubjectFolder, 'Results', 'Profiles', 'Surfaces');
    
    for ihs = 1:2
        
        %% ROI
%         cd(fullfile(SubjectFolder, 'Structural', 'CBS'))
%         vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' ...
%             hs(ihs) 'cr_gm_avg_inf.vtk']);
%         [vertex,faces,~] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', vtk.name), 0, 1);
%         
%         
%         cd(fullfile(SubjectFolder, 'ROI_MIPAV'));
%         [~, ~, mapping] = read_vtk(['Subj_' SubjID '_' hs(ihs) 'cr_AllROIs.vtk'], 0, 1);
%         write_vtk(['Subj_' SubjID '_'  hs(ihs) 'cr_AllROIs_inf.vtk'], vertex, faces, mapping)
%         
%         g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_AllROIs_inf.vtk']);
%         
%         g.cdata(g.cdata==10)=3;
%         g.cdata(g.cdata==15)=4;
%         
%         close all
%         Fig = figure('name', 'ROI' ,'Visible', 'off');
%         colormap([0.8 0.8 0.8;1 0 0; 0 0 1;1 0 0; 0 0 1])
%         plot(g,g);
%         
%         PrintFigSurf(SubjID, Fig, hs, ihs, 'ROI_inf', '', Sufix, [0 4])
        
        
        for iVar=1:numel(Var)
            
            
            %% Sensory modalities
            for iSens=1:numel(Sens)
                
                if Smooth
                    cd(Results_Folder)
                else
                    cd(fullfile(Results_Folder,'Cdtions'))
                end
                
                if Smooth
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' Sens{iSens} '_' Var{iVar} '_smoothdata.vtk']);
                else
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' Sens{iSens} '_' Var{iVar} '.vtk']);
                end
                
                close all
                Fig = figure('name', '1' ,'Visible', 'off');
                colormap([R' G' B']);
                
                plot(g,g);
                
                PrintFigSurf(SubjID, Fig, hs, ihs, Sens{iSens}, Var{iVar}, Sufix, CLIM{iVar})
                
            end
            
            
            %% CrossSensory modalities
            for iCrossSens=1:numel(CrossSens)
                
                if Smooth
                    cd(Results_Folder)
                else
                    cd(fullfile(Results_Folder,'CrossSens'))
                end
                
                if Smooth
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' CrossSens{iCrossSens} '_' Var{iVar} '_smoothdata.vtk']);
                else
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' CrossSens{iCrossSens} '_' Var{iVar} '.vtk']);
                end
                
                close all
                Fig = figure('name', '1' ,'Visible', 'off');
                colormap([R' G' B']);
                
                plot(g,g);
                
                PrintFigSurf(SubjID, Fig, hs, ihs, CrossSens{iCrossSens}, Var{iVar}, Sufix, CLIM{iVar})
                
                
            end
            
            
            %% Attention
            for iAtt=1:numel(Att)
                
                if Smooth
                    cd(Results_Folder)
                else
                    cd(fullfile(Results_Folder,'Att'))
                end
                
                if Smooth
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' Att{iAtt} '_' Var{iVar} '_smoothdata.vtk']);
                else
                    g = gifti(['Subj_' SubjID '_' hs(ihs) 'cr_' Att{iAtt} '_' Var{iVar} '.vtk']);
                end
                
                close all
                Fig = figure('name', '1' ,'Visible', 'off');
                colormap([R' G' B']);
                
                plot(g,g);
                
                PrintFigSurf(SubjID, Fig, hs, ihs, Att{iAtt}, Var{iVar}, Sufix, CLIM{iVar})
                
            end
            
        end
    end
    
end

