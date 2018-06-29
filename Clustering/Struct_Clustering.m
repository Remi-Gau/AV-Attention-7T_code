
%%
clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

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


FontSize = 10;

NbLayers = 10;

MaxNbCluster = 3;

Replicates = 200;

MaxIter = 200;

Visibility = 'on';

ClusterColors = [...
    '.b';
    '.r';
    '.g';
    '.c';
    '.m';
    '.k';
    '.y'];


for SubjInd = 5 %:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    cd(fullfile(StartFolder, 'Subjects_Data'))
    
    for hs=1:2
        
        if hs==1
            fprintf(' Running clustering on left HS\n')
            Suffix = 'lcr';
        else
            fprintf(' Running clustering on right HS\n')
            Suffix = 'rcr';
        end
        
        %% Get data
        
%         cd(fullfile(StartFolder, 'Subjects_Data', 'MNI_ROI'))    
%         ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_L.vtk']);            
%         [~,~,mapping] = read_vtk(ROI_vtk.name, 0, 1);       
%         VertOfInt = find(any([mapping==101;mapping==102;mapping==103]));
        
        cd(fullfile(StartFolder, 'Subjects_Data'))
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te10.nii.vtk']);
        [~,~,TE10] = read_vtk(ROI_vtk.name, 0, 1);
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te11.nii.vtk']);
        [~,~,TE11] = read_vtk(ROI_vtk.name, 0, 1);
        ROI_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_wmAuditory_Te12.nii.vtk']);
        [~,~,TE12] = read_vtk(ROI_vtk.name, 0, 1);
        
        VertOfInt = find(any([TE10;TE11;TE12]>.1));       
    
        cd(fullfile(StartFolder, 'Subjects_Data', 'Clustering'))
        if hs==1
        qT1Map_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_T1_10_Layers.vtk']);
        else
            qT1Map_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_T1_10_Layers_R.vtk']);
        end

        [~,~,Profiles] = read_vtk(qT1Map_vtk.name, NbLayers, 1);
        
        Profiles = Profiles(:,VertOfInt)';
        
        
        
        %% Clustering        
        inf_vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' Suffix '_gm_avg_inf.vtk']);       

        [vertex,face,mapping] = read_vtk(inf_vtk.name, 0, 1);  
        
        idx = [];
        C = {};
        sumd = [];
        
        for NbCluster=3 %:MaxNbCluster
            fprintf('  Running clustering with %i clusters\n', NbCluster)
            [idx(:,NbCluster),C{NbCluster},sumd] = kmeans(Profiles,NbCluster, ...
                'EmptyAction', 'drop', 'Replicates', Replicates, 'MaxIter', MaxIter); %#ok<SAGROW>
            
            mapping = zeros(size(mapping));
            mapping(VertOfInt) = idx(:,NbCluster);
            
            if NbCluster>1
                write_vtk(['T1_' SubjID '_' Suffix '_Cluster_' num2str(NbCluster) '.vtk'], vertex, face, mapping)
            end

            TotalDist(NbCluster) = sum(sumd);

        end
        
        
        
        %%
%         H = figure('Name', ['Scree plot: ' num2str(NbLayers) ' layers'], 'Position', [500, 500, 750, 750], ...
%             'Color', [1 1 1], 'Visible', Visibility);
%         plot(TotalDist, 'linewidth', 2)
%         
%         t=ylabel('Total distance');
%         set(t,'fontsize', FontSize);
%         
%         t=xlabel('Number of clusters');
%         set(t,'fontsize', FontSize);
%         
%         set(gca,'tickdir', 'out', ...
%             'xtick', 1:length(TotalDist), ...
%             'xticklabel', 1:MaxNbCluster, ...
%             'ticklength', [0.01 0.01], 'fontsize', FontSize)
%         
%         if hs==1
%             print(gcf, ['Scree_plot_lcr_Subj' SubjID '_' num2str(NbLayers) '_Layers.tif'], '-dtiff')
%         else
%             print(gcf, ['Scree_plot_rcr_Subj' SubjID '_' num2str(NbLayers) '_Layers.tif'], '-dtiff')
%         end
%         
%         
%         NbCluster = find(diff(diff(TotalDist))==max( diff(diff(TotalDist))))+2;
%         
%         idx = idx(:,NbCluster);
%         C = C{NbCluster};
%         
%         clear t H TotalDist



        %%
%         figure('Name', 'Clusters', 'Position', [500, 500, 1000, 1000], ...
%             'Color', [1 1 1], 'Visible', Visibility);
%         
%         COMB = combnk(1:size(C,2),2);
%         SubplotTarget = max(COMB(:)-1)*(COMB(:,1)-1)+COMB(:,2)-1;
%         
%         for iLayer=1:size(COMB,1)
%             
%             subplot(size(C,2)-1,size(C,2)-1,SubplotTarget(iLayer))
%             hold on
%             for iCluster=1:NbCluster
%                 plot(Profiles(idx==iCluster, COMB(iLayer,1)), ...
%                     Profiles(idx==iCluster, COMB(iLayer,2)),...
%                     ClusterColors(iCluster,:))
%             end
%             t=xlabel(['Layer ' num2str(COMB(iLayer,1))]);
%             set(t,'fontsize', FontSize);
%             
%             t=ylabel(['Layer ' num2str(COMB(iLayer,2))]);
%             set(t,'fontsize', FontSize);
%             
%                 set(gca,'tickdir', 'out', ...
%                     'ytick', [0 10000], ...
%                     'yticklabel', [0 10000], ...
%                     'xtick', [0 10000], ...
%                     'xticklabel', [0 10000], ...
%                     'ticklength', [0.01 0.01], 'fontsize', FontSize)
%         end
%         
%         if hs==1
%             print(gcf, ['Clusters_lcr_Subj_' SubjID '_' num2str(NbLayers) '_Layers.tif'], '-dtiff')
%         else
%             print(gcf, ['Clusters_rcr_Subj_' SubjID '_' num2str(NbLayers) '_Layers.tif'], '-dtiff')
%         end
%         
%         clear t COMB iLayer
        


        %%
%         figure('Name', 'Clusters profiles', 'Position', [500, 500, 1000, 1000], ...
%             'Color', [1 1 1], 'Visible', 'on');
%         
%         for iCluster=1:NbCluster
%             subplot(1,NbCluster,iCluster)
%             plot(C(iCluster,:))
%         end
%         
        
        
    end
    
    
end


cd(StartFolder)