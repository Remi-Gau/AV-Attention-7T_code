clear; clc;

more off

% NbLayers = 4;
% WithQuad = 1;
% WithPerm = 1;

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




for NbLayers=[6]
    for WithQuad= [0]
        for WithPerm = [0]
   
            
            FigureFolder = fullfile(StartFolder, 'Figures', 'Profiles');
            cd(FigureFolder)
            if WithQuad
                load(fullfile(FigureFolder, [num2str(NbLayers) '_layers'], strcat('Data_Block_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data') %#ok<*UNRCH>
            else
                load(fullfile(FigureFolder, [num2str(NbLayers) '_layers'], strcat('Data_Block_NoQuadGLM_', num2str(NbLayers), '_Layers', '.mat')), 'AllSubjects_Data')
            end
            
            
            %%
            for iCond = 1:2
                for iROI=1:length(AllSubjects_Data)
                    Include = find(AllSubjects_Data(iROI).Include);
                    tmp = squeeze(AllSubjects_Data(iROI).Differential.Blocks.MainEffects.Beta.DATA(:,iCond,Include))';
                    [~,P] = ttest(tmp);
                    All_P(iROI,:,iCond) = P;
                    
                    All_Betas(:,iROI,1:size(tmp,2),iCond) = tmp;
                end
            end
            clear P
            
            
            %%
            HolmsThresholds = .05*ones(1,length(AllSubjects_Data))./(length(AllSubjects_Data)+1-(1:length(AllSubjects_Data)));
            
            for iCond = 1:2
                for iP=1:size(All_P,2)
                    tmp = All_P(:,iP,iCond);
                    [~,I] = ismember(tmp,sort(tmp,'ascend'));
                    Thresholds(:,iP,iCond) = HolmsThresholds(I); %#ok<*SAGROW>
                end
            end
            clear All_P
            
            
            %% Permutations
            % Get all possible permutations
            for iSubj=1:size(All_Betas,1)
                sets{iSubj} = [0 1];
            end
            [a,b,c,d,e,f,g,h,i,j,k] = ndgrid(sets{:});
            cartProd = logical([a(:),b(:),c(:),d(:),e(:),f(:),g(:),h(:),i(:),j(:),k(:)]);
            
            % Do the permutations
            NbPerm=size(cartProd,1);
            NullDist = nan(NbPerm,size(All_Betas,3));
            
            for iCond = 1:2
                for iPerm=1:NbPerm
                    tmp = All_Betas(:,:,:,iCond);
                    tmp(cartProd(iPerm,:),:,:) = tmp(cartProd(iPerm,:),:,:)*-1;
                    NullDist(iPerm,:) = squeeze(max(abs(mean(tmp)),[],2));
                end
                AllNullDist(:,:,iCond) = sort(NullDist);
                clear NullDist
            end
            clear All_Betas
            
            
            
            %% Plots
            for ROI_Ind = 1:length(AllSubjects_Data)
                
                close all
                
                fprintf([AllSubjects_Data(ROI_Ind).name '\n'])
                
                Include = find(AllSubjects_Data(ROI_Ind).Include);
                
                if ~isempty(Include)
                    
                    %% Plot All +/- STD with subjects
                    Visible='on';
                    PRINT=2;
                    
                    if WithPerm
                        Name = [AllSubjects_Data(ROI_Ind).name '-Perm'];
                    else
                        Name = [AllSubjects_Data(ROI_Ind).name '-NoPerm'];
                    end
                    
                    if WithQuad
                        Name = [Name '-ALL'];
                    else
                        Name = [Name '-NoQuad-ALL'];
                    end
                    
                    ToPlot.shadedErrorBar.mean = cat(2,flipud(AllSubjects_Data(ROI_Ind).MEAN), ...
                        flipud(AllSubjects_Data(ROI_Ind).Differential.MainEffect.MEAN));
                    
                    ToPlot.shadedErrorBar.errorbar = cat(2,flipud(AllSubjects_Data(ROI_Ind).SEM),...
                        flipud(AllSubjects_Data(ROI_Ind).Differential.MainEffect.SEM));
                    
                    %         ToPlot.shadedErrorBar.errorbar = cat(2,flipud(AllSubjects_Data(ROI_Ind).STD),...
                    %             flipud(AllSubjects_Data(ROI_Ind).Differential.MainEffect.STD));
                    
                    ToPlot.Subjects.Data = AllSubjects_Data(ROI_Ind).DATA(:,:,Include);
                    ToPlot.Max(1) = max(ToPlot.Subjects.Data(:));
                    ToPlot.Min(1) = min(ToPlot.Subjects.Data(:));
                    
                    tmp = AllSubjects_Data(ROI_Ind).Differential.MainEffect.DATA;
                    ToPlot.Max(2:3) = max(max(tmp,[],3));
                    ToPlot.Min(2:3) = min(min(tmp,[],3));
                    
                    ToPlot.Subjects.Data = cat(2,ToPlot.Subjects.Data, tmp);
                    
                    ToPlot.Include=Include;
                    
                    ToPlot.WithQuad = WithQuad;
                    
                    if WithPerm
                        ToPlot.AllNullDist=AllNullDist;
                    else
                        ToPlot.Thresholds=squeeze(Thresholds(ROI_Ind,:,:));
                    end
                    
                    ToPlot.Beta = cat(2,AllSubjects_Data(ROI_Ind).Differential.Blocks.Beta.DATA(:,[1:3 5:7],Include),...
                        AllSubjects_Data(ROI_Ind).Differential.Blocks.MainEffects.Beta.DATA(:,:,Include));
                    
                    PlotLayersForFig(ToPlot, Name, Visible)
                    
                    clear ToPlot
                    
                    
                end
                
            end
            cd(StartFolder)
            clear AllNullDist
        end
    end
end


