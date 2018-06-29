%%
clear; clc; close all;

StartFolder = fullfile(pwd, '..', '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')))

Visible = 'off';

Fontsize = 8;

NbLayers = 10;

Bins = 500;

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
    '16'
    ];

mn = size(SubjectList,1);
n  = round(mn^0.4);
m  = ceil(mn/n);


%%
for SubjInd = 1%:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);
    
    cd(fullfile(StartFolder, 'Subjects_Data'))
    
    %%
    for hs=1%:2
        
        if hs==1
            fprintf(' Doing left HS\n')
            Suffix = 'lcr';
        else
            fprintf(' Doing right HS\n')
            Suffix = 'rcr';
        end
        
        
        %% Get data
        cd(fullfile(SubjectFolder, 'Structural', 'CBS'))
        inf_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf.vtk']);
        [vertex,face,mapping] = read_vtk(inf_vtk.name, 0, 1);
        
        cd(fullfile(SubjectFolder, 'Structural', 'CBS'))
        qT1Map_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_qT1.vtk']);
        [~,~,qT1MapOri] = read_vtk(qT1Map_vtk.name, 0, 1);
        
        
        cd(fullfile(StartFolder, 'Subjects_Data', 'Curvature'))
        if hs==1
            qT1Map10_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_T1_10_Layers.vtk']);
        else
            qT1Map10_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_T1_10_Layers_R.vtk']);
        end
        [~,~,qT1Map10] = read_vtk(qT1Map10_vtk.name, 0, 1);
        
        MeanCurvature_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_MeanCurvature.vtk']);
        [~,~,MeanCurvature] = read_vtk(MeanCurvature_vtk.name, 0, 1);
        
        GaussCurvature_vtk = dir(['T1_' SubjID '*' Suffix '_gm_avg_inf_GaussCurvature.vtk']);
        [~,~,GaussCurvature] = read_vtk(GaussCurvature_vtk.name, 0, 1);
        
        
        %%
        figure(1);
        clf
        
%         Yid=find(MeanCurvature~=0);
        Yid =1:numel(MeanCurvature);
        [B,DEV,STATS] = glmfit([MeanCurvature(Yid)' ones(size(Yid'))] , qT1MapOri(Yid), 'normal', 'constant', 'off');
        mapping = zeros(size(MeanCurvature));
        mapping(Yid) = B(2)+STATS.resid;
        write_vtk(['T1_' SubjID '_' Suffix  '_qT1_MeanCurvatureDetrend.vtk'], vertex, face, mapping)
        
        subplot(221)
        hold on
        scatter(MeanCurvature(Yid), qT1MapOri(Yid), '.')
        plot([min(MeanCurvature(Yid)) max(MeanCurvature(Yid))], [min(MeanCurvature(Yid)) max(MeanCurvature(Yid))]*B(1)+B(2), 'r')
        title('MeanCurvature')
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('MeanCurvature (mm)');
        set(t,'fontsize',Fontsize);
        axis([-1 1 0 4000])
        
        subplot(223)
        hold on
        Xi = linspace(1,500,Bins);
        Yi = linspace(1,4000,Bins);
        
        Xr = interp1(Xi, 1:numel(Xi), (MeanCurvature(Yid)+2.5)*Bins/5, 'nearest');
        Yr = interp1(Yi, 1:numel(Yi), qT1MapOri(Yid), 'nearest');
        
        Yr(isnan(Xr))=[];
        Xr(isnan(Xr))=[];
        
        Xr(isnan(Yr))=[];
        Yr(isnan(Yr))=[];
        
        PDF = accumarray([round(Yr') round(Xr')], 1, [Bins Bins]);
        PDF=PDF/max(PDF(:));
        
        B = glmfit(Xr, Yr, 'normal');
        
        
        A=colormap;
        A(1,:) = [1 1 1];
        colormap(A);
        imagesc(PDF)
        axis xy
        
        plot([0 500], ([0 500])*B(2)+B(1), 'k', 'LineWidth', 1)
        
        axis([175 326 150 450])
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
            'xtick', linspace(1, 500, 21) ,'xticklabel', linspace(-2.5, 2.5, 21),...
            'ytick', linspace(1, 500, 6) ,'yticklabel', linspace(0, 4000, 6))
        
        
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('MeanCurvature (mm)');
        set(t,'fontsize',Fontsize);
        
        
        
        
        
%         Yid=find(all([GaussCurvature~=0;GaussCurvature<2.5;GaussCurvature>-2.5]));
         Yid =1:numel(GaussCurvature);
        [B,DEV,STATS] = glmfit([GaussCurvature(Yid)' ones(size(Yid'))] , qT1MapOri(Yid), 'normal', 'constant', 'off');
        mapping = zeros(size(MeanCurvature));
        mapping(Yid) = B(2)+STATS.resid;
        write_vtk(['T1_' SubjID '_' Suffix  '_qT1_GaussCurvatureDetrend.vtk'], vertex, face, mapping)
        
        subplot(222)
        hold on
        scatter(GaussCurvature(Yid), qT1MapOri(Yid), '.')
        plot([min(GaussCurvature(Yid)) max(GaussCurvature(Yid))], [min(GaussCurvature(Yid)) max(GaussCurvature(Yid))]*B(1)+B(2), 'r')
        title('GaussCurvature')
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('GaussCurvature (mm)');
        set(t,'fontsize',Fontsize);%
        axis([-1 1 0 4000])
        
        subplot(224)
        NbBins = 1000;
        hold on
        Xi = linspace(1,1000,NbBins);
        Yi = linspace(1,4000,NbBins);
        
        Xr = interp1(Xi, 1:numel(Xi), (GaussCurvature(Yid)+2.5)*NbBins/5, 'nearest');
        Yr = interp1(Yi, 1:numel(Yi), qT1MapOri(Yid), 'nearest');
        
        Yr(isnan(Xr))=[];
        Xr(isnan(Xr))=[];
        
        Xr(isnan(Yr))=[];
        Yr(isnan(Yr))=[];
        
        PDF = accumarray([round(Yr') round(Xr')], 1, [NbBins NbBins]);
        PDF=PDF/max(PDF(:));
        
        B = glmfit(Xr, Yr, 'normal');
        
        
        A=colormap;
        A(1,:) = [1 1 1];
        colormap(A);
        imagesc(PDF)
        axis xy
        
        plot([0 NbBins], ([0 NbBins])*B(2)+B(1), 'k', 'LineWidth', 1)
        
        axis([450 550 300 800])
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
            'xtick', linspace(1, NbBins, 41) ,'xticklabel', linspace(-2.5, 2.5, 41),...
            'ytick', linspace(1, NbBins, 6) ,'yticklabel', linspace(0, 4000, 6))
        
        
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('GaussCurvature (mm)');
        set(t,'fontsize',Fontsize);%
        
        print(gcf, ['Curvature_Subj' SubjID '_' Suffix '.tif'] ,'-dtiff')
        
        
    end
end
