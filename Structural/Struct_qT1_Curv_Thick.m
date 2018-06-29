%%
clc; clear; close all;

StartDirectory = pwd;

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

XP = [...
    'exp-0000';...
    'exp-0001';...
    'exp-0002';...
    'exp-0003';...
    'exp-0004';...
    'exp-0005';...
    'exp-0006';...
    'exp-0007';...
    'exp-0008';...
    'exp-0009';...
    'exp-0010';...
    'exp-0011';...
    'exp-0012'
    ];

Visible = 'off';

set(0, 'Visible', Visible)

Fontsize = 8;

Bins = 500;

mn = size(SubjectList,1);
n  = round(mn^0.4);
m  = ceil(mn/n);


%%
for SubjInd = 1:size(SubjectList,1)
    
    SubjID = SubjectList(SubjInd,:);
    
    SubjectFolder = fullfile(StartDirectory, 'Subjects_Data', ['Subject_' SubjID]);
    
    SubXP = XP(SubjInd,:);
    
    fprintf(['\nProcessing subject : ' SubjID '\n'])
    
    
    % LEFT
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'ThickCurv', SubXP, [SubXP '-AA'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*lcr*inf_thickness.vtk']));
    [~, ~, Thick_L] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    clear LogFileList;
    
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'ThickCurv', SubXP, [SubXP '-BA'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*lcr*inf_mean.vtk']));
    [~, ~, Curv_L] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    clear LogFileList;
    
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'Vertex', 'T1', ...
        SubXP, [SubXP '-CB'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*lcr*inf*10*.vtk']));
    [~, ~, T1_mapping_L] = read_vtk(fullfile (pwd, LogFileList.name), 10, 1);
    clear LogFileList;
    
    
    % RIGHT
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'ThickCurv', SubXP, [SubXP '-AB'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*rcr*inf_thickness.vtk']));
    [~, ~, Thick_R] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    clear LogFileList;
    
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'ThickCurv', SubXP, [SubXP '-BB'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*rcr*inf_mean.vtk']));
    [~, ~, Curv_R] = read_vtk(fullfile (pwd, LogFileList.name), 0, 1);
    clear LogFileList;
    
    cd(fullfile(SubjectFolder, 'Structural', 'CBS', 'Vertex', 'T1', ...
        SubXP, [SubXP '-CA'], 'SurfaceMeshMapping'))
    LogFileList = dir(strcat(['T1_' SubjID '*rcr*inf*10*.vtk']));
    [~, ~, T1_mapping_R] = read_vtk(fullfile (pwd, LogFileList.name), 10, 1);
    clear LogFileList;
    
    
    for DepthInd=1:size(T1_mapping_L,1)
        
        T1 = [T1_mapping_L(DepthInd,:) T1_mapping_R(DepthInd,:)];
        Thick = [Thick_L Thick_R];
        
        Thick(any([T1<10;T1>3900])) = [];
        T1(any([T1<10;T1>3900])) = [];
        T1(Thick==0) = [];
        Thick(Thick==0) = [];
        
        [r,p]=corrcoef(Thick,T1);
        R.Thick(DepthInd,SubjInd)=r(2,1).^2;
        P.Thick(DepthInd,SubjInd) = p(2,1); clear r p
        
        
        
        %% Plot Thick VS T1
        Xi = linspace(1,500,Bins);
        Yi = linspace(1,4000,Bins);
        
        Xr = interp1(Xi, 1:numel(Xi), Thick*Bins/5, 'nearest'); clear Xi
        Yr = interp1(Yi, 1:numel(Yi), T1, 'nearest'); clear Yi
        
        Yr(isnan(Xr))=[];
        Xr(isnan(Xr))=[];
        
        Xr(isnan(Yr))=[];
        Yr(isnan(Yr))=[];
        
        PDF = accumarray([round(Yr) ; round(Xr)]', 1, [Bins Bins]);
        PDF=PDF/max(PDF(:));
        PDF_all{DepthInd,SubjInd,1} = PDF;
        
        B = glmfit(Xr, Yr, 'normal'); clear Yr Xr
        
        if SubjInd==1
            H_Thick(DepthInd) = figure('Name', ['Thickness VS T1 at depth: ' sprintf('%3.1f', DepthInd)], ...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
        elseif SubjInd~=1
            figure(H_Thick(DepthInd))
            set(H_Thick(DepthInd),'Visible', Visible)
        end
        
        subplot(m,n,SubjInd)
        hold on
        
        A=colormap;
        A(1,:) = [1 1 1];
        colormap(A);
        imagesc(PDF); clear PDF
        axis xy
        
        plot([0 500], ([0 500])*B(2)+B(1), 'k', 'LineWidth', 1); clear B
        
        axis([0 500 100 500])
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
            'xtick', linspace(1, 500, 9) ,'xticklabel', linspace(0, 6, 9),...
            'ytick', linspace(1, 500, 6) ,'yticklabel', linspace(0, 4000, 6))
        
        t=title(['Subject ' SubjectList(SubjInd,:)]);
        set(t,'fontsize', Fontsize);
        
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('Thickness (mm)');
        set(t,'fontsize',Fontsize); clear t
        
        
        %% Plot Curv VS T1
        T1 = [T1_mapping_L(DepthInd,:) T1_mapping_R(DepthInd,:)];
        Curv = [Curv_L Curv_R];
        
        Curv(any([T1<10;T1>3900])) = [];
        T1(any([T1<10;T1>3900])) = [];
%         T1(Curv>1) = [];
%         Curv(Curv>1) = [];
        
        [r,p]=corrcoef(Curv,T1);
        R.Curv(DepthInd,SubjInd)=r(2,1).^2;
        P.Curv(DepthInd,SubjInd) = p(2,1); clear r p
        
        Xi = linspace(1,500,Bins);
        Yi = linspace(1,4000,Bins);
        
        Xr = interp1(Xi, 1:numel(Xi), (Curv)*Bins/5, 'nearest'); clear Curv_tmp
        Yr = interp1(Yi, 1:numel(Yi), T1, 'nearest');
        
        Yr(isnan(Xr))=[];
        Xr(isnan(Xr))=[];
        
        Xr(isnan(Yr))=[];
        Yr(isnan(Yr))=[];
        
        PDF = accumarray([round(Yr) ; round(Xr)]', 1, [Bins Bins]);
        PDF=PDF/max(PDF(:));
        PDF_all{DepthInd,SubjInd,2} = PDF;
        
        B = glmfit(Xr, Yr, 'normal');
        
        if SubjInd==1
            H_Curv(DepthInd) = figure('Name', ['Curvature VS T1 at depth: ' sprintf('%3.1f', DepthInd)], ...
                'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', Visible);
        elseif SubjInd~=1
            figure(H_Curv(DepthInd));
            set(H_Curv(DepthInd),'Visible', Visible)
        end
        
        subplot(m,n,SubjInd)
        hold on
        
        A=colormap;
        A(1,:) = [1 1 1];
        colormap(A);
        imagesc(PDF)
        axis xy
        
        plot([0 500], ([0 500])*B(2)+B(1), 'k', 'LineWidth', 2)
        
        axis([0 75 150 450])
        
        set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
            'xtick', linspace(1, 500, 21) ,'xticklabel', linspace(0, 1, 21),...
            'ytick', linspace(1, 500, 6) ,'yticklabel', linspace(0, 4000, 6))
        
        t=title(['Subject ' SubjectList(SubjInd,:)]);
        set(t,'fontsize', Fontsize);
        
        t=ylabel('T1 value');
        set(t,'fontsize',Fontsize);
        
        t=xlabel('Curvature (Deg)');
        set(t,'fontsize',Fontsize);
        
    end
    
    clear Curv Thick T1 T1_mapping_R Thick_R Curv_R T1_mapping_L Thick_L Curv_L
    
end

cd(StartDirectory)

save('CorrT1CurvThick.mat')


%% Print figures
cd('Vertices_Analysis')
for DepthInd = 1:length(H_Curv)
    print(H_Curv(DepthInd), ['Curvature VS T1 at depth: ' sprintf('%3.1f', DepthInd) '.tif'], '-dtiff')
    print(H_Thick(DepthInd), ['Thickness VS T1 at depth: ' sprintf('%3.1f', DepthInd) '.tif'], '-dtiff')
end

figure(H_Curv(6))

figure(H_Thick(6))


%%
figure('Name', 'R^2 values', ...
    'Position', [100, 100, 1500, 1000], 'Color', [1 1 1], 'Visible', 'on');
subplot(211)
errorbar(1:size(R.Curv,1), mean(R.Curv,2), nansem(R.Curv,2), 'k', 'LineWidth', 2)

set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
    'xtick', 1:11 ,'xticklabel', 0:0.1:1)

t=ylabel('R^2 between curvature and T1');
set(t,'fontsize',Fontsize);

t=xlabel('Depth fraction');
set(t,'fontsize',Fontsize);

subplot(212)
errorbar(1:size(R.Thick,1), mean(R.Thick,2), nansem(R.Thick,2), 'k', 'LineWidth', 2)

t=ylabel('R^2 between thickness and T1');
set(t,'fontsize',Fontsize);

t=xlabel('Depth fraction');

set(t,'fontsize',Fontsize);

set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', Fontsize, ...
    'xtick', 1:size(R.Thick,1) ,'xticklabel', linspace(0,1,size(R.Thick,1)))

print(gcf, 'R2 as function of depth.tif', '-dtiff')