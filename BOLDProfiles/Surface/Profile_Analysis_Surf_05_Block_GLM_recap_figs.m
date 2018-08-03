%%
close all
clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

FigureFolder = fullfile(StartFolder,'Figures','Profiles','Surfaces');
mkdir(FigureFolder)

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
CLIM = {[-2 2], [-1 1]};
Var={'Cst','Lin'};
ROI = {'_A1-PT','_V1','_V2-3'};
ROI_folder = {'A1','V1','V2-3'};
Smooth = 1;

if Smooth
    Sufix = '_smooth_6fwhm';
    SmoothFolder = 'Smooth6';
else
    Sufix = ''; %#ok<*UNRCH>
    SmoothFolder = 'NonSmooth';
end

FigDim = [100, 100, 1500, 1500];
Visibility = 'off';


%% Basic
fprintf('Basic\n')

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    fprintf(' Subject %s\n',SubjID)
    
    close all
    
    for ihs = 1:2
        
        for iROI=1:numel(ROI)
            
            Fig=figure('Name', 'Basic conditions', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            
            set(Fig, 'PaperUnits','centimeters');
            set(Fig, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(Fig, 'PaperPositionMode', 'manual');
            set(Fig, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            
            iSubplot = 1;
            
            for iSens=1:numel(Sens)
                
                for iVar=1:numel(Var)
                    subplot(3,3,iSubplot)
                    
                    SubjID = SubjectList(SubjInd,:);
                    if Smooth
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces');
                    else
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces','Cdtions');
                    end
                    Img = imread(fullfile(Folder, ['Subj_' SubjID '_' hs(ihs) 'cr_' Sens{iSens} '_' Var{iVar} ROI{iROI} Sufix '.tif']));
                    imagesc(Img)
                    %         axis off
                    set(gca, 'Xtick', [], 'Ytick', [])
                    axis image
                    
                    iSubplot = iSubplot + 1;
                end
                
                subplot(3,3,iSubplot)
                Img = imread(fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
                    'ROI_MIPAV' ,['Subj_' SubjID '_' hs(ihs) 'cr_ROI_inf_' ROI{iROI} '.tif']));
                imagesc(Img)
                set(gca, 'Xtick', [], 'Ytick', [])
                axis image
                
                iSubplot = iSubplot + 1;
                
                
            end
            
            subplot(3,3,1)
            title('Constant')
            t = ylabel('Audio');
            set(t,'fontsize',12)
            
            subplot(3,3,2)
            title('Linear')
            
            subplot(3,3,3)
            title('ROI')
            
            subplot(3,3,4)
            t = ylabel('Visual');
            set(t,'fontsize',12)
            
            subplot(3,3,7)
            t = ylabel('Audio-visual');
            set(t,'fontsize',12)
            
            mtit(['Subject ' SubjID ' - ' hs(ihs) 'cr - ' ROI{iROI}(2:end)], 'fontsize', 14, 'xoff',0,'yoff',.025)
            set(Fig,'Visible',Visibility)
            
            print(Fig, fullfile(FigureFolder,'Baseline', ROI_folder{iROI}, SmoothFolder,...
                ['Subj_' SubjID '_' hs(ihs) 'cr_Baseline_'  ROI{iROI} Sufix '.tif']), '-dtiff')
            
        end
        
    end
    
end



%% CrossSens
fprintf('CrossSens\n')

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    fprintf(' Subject %s\n',SubjID)
    
    close all
    
    for ihs = 1:2
        
        for iROI=1:numel(ROI)
            
            Fig=figure('Name', 'CrossSens', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            
            set(Fig, 'PaperUnits','centimeters');
            set(Fig, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(Fig, 'PaperPositionMode', 'manual');
            set(Fig, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            
            iSubplot = 1;
            
            for iSens=1:numel(CrossSens)
                
                for iVar=1:numel(Var)
                    subplot(2,3,iSubplot)
                    
                    SubjID = SubjectList(SubjInd,:);
                    if Smooth
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces');
                    else
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces','CrossSens');
                    end
                    Img = imread(fullfile(Folder, ['Subj_' SubjID '_' hs(ihs) 'cr_' CrossSens{iSens} '_' Var{iVar} ROI{iROI} Sufix '.tif']));
                    imagesc(Img)
                    
                    set(gca, 'Xtick', [], 'Ytick', [])
                    axis image
                    
                    iSubplot = iSubplot + 1;
                end
                
                subplot(2,3,iSubplot)
                Img = imread(fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
                    'ROI_MIPAV' ,['Subj_' SubjID '_' hs(ihs) 'cr_ROI_inf_' ROI{iROI} '.tif']));
                imagesc(Img)
                set(gca, 'Xtick', [], 'Ytick', [])
                axis image
                
                iSubplot = iSubplot + 1;
                
            end
            
            subplot(2,3,1)
            title('Constant')
            t = ylabel('AV-A');
            set(t,'fontsize',12)
            
            subplot(2,3,2)
            title('Linear')
            
            subplot(2,3,3)
            title('ROI')
            
            subplot(2,3,4)
            t = ylabel('AV-V');
            set(t,'fontsize',12)
            
            mtit(['Subject ' SubjID ' - ' hs(ihs) 'cr - ' ROI{iROI}(2:end)], 'fontsize', 14, 'xoff',0,'yoff',.025)
            set(Fig,'Visible',Visibility)
            
            print(Fig, fullfile(FigureFolder,'CrossSensory', ROI_folder{iROI}, SmoothFolder, ...
                ['Subj_' SubjID '_' hs(ihs) 'cr_CrossSensory_'  ROI{iROI} Sufix '.tif']), '-dtiff')
            
        end
        
    end
    
end



%% Attention
fprintf('Attention\n')

for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    fprintf(' Subject %s\n',SubjID)
    
    close all
    
    for ihs = 1:2
        
        for iROI=1:numel(ROI)
            
            Fig=figure('Name', 'Attention', 'Position', FigDim, 'Color', [1 1 1], 'Visible', Visibility);
            
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            
            set(Fig, 'PaperUnits','centimeters');
            set(Fig, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(Fig, 'PaperPositionMode', 'manual');
            set(Fig, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            
            iSubplot = 1;
            
            for iAtt=1:numel(Att)
                
                for iVar=1:numel(Var)
                    subplot(4,3,iSubplot)
                    
                    SubjID = SubjectList(SubjInd,:);
                    if Smooth
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces');
                    else
                        Folder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], 'Results', 'Profiles', 'Surfaces','Att');
                    end
                    Img = imread(fullfile(Folder, ['Subj_' SubjID '_' hs(ihs) 'cr_' Att{iAtt} '_' Var{iVar} ROI{iROI} Sufix '.tif']));
                    imagesc(Img)
                    
                    set(gca, 'Xtick', [], 'Ytick', [])
                    axis image
                    
                    iSubplot = iSubplot + 1;
                end
                
                subplot(4,3,iSubplot)
                Img = imread(fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID], ...
                    'ROI_MIPAV' ,['Subj_' SubjID '_' hs(ihs) 'cr_ROI_inf_' ROI{iROI} '.tif']));
                imagesc(Img)
                set(gca, 'Xtick', [], 'Ytick', [])
                axis image
                
                iSubplot = iSubplot + 1;
                
            end
            
            subplot(4,3,1)
            title('Constant')
            t = ylabel('A Stim: A Att-V Att');
            set(t,'fontsize',12)
            
            subplot(4,3,2)
            title('Linear')
            
            subplot(4,3,3)
            title('ROI')
            
            subplot(4,3,4)
            t = ylabel('V Stim: A Att-V Att');
            set(t,'fontsize',12)
            
            subplot(4,3,7)
            t = ylabel('AV Stim: A Att-V Att');
            set(t,'fontsize',12)
            
            subplot(4,3,10)
            t = ylabel('A Att-V Att');
            set(t,'fontsize',12)
            
            mtit(['Subject ' SubjID ' - ' hs(ihs) 'cr - ' ROI{iROI}(2:end)], 'fontsize', 14, 'xoff',0,'yoff',.025)
            set(Fig,'Visible',Visibility)
            
            print(Fig, fullfile(FigureFolder,'Attention', ROI_folder{iROI}, SmoothFolder,...
                ['Subj_' SubjID '_' hs(ihs) 'cr_Attention_'  ROI{iROI} Sufix '.tif']), '-dtiff')
            
        end
        
    end
    
end
