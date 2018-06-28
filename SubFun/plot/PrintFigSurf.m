function PrintFigSurf(SubjID, Fig, hs, ihs, Modality, Var, Sufix, CLIM)
axis off

axis1 = gca;

set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1.66 1.31 1],'Projection','perspective', 'CLIM', CLIM);

%%
if ihs==1
    set(gca,'CameraPosition',[2100 1800 -200],'CameraTarget',[140 115 89.5],...
        'CameraUpVector',[-0.08 -0.07 -0.99],'CameraViewAngle',1);
    axis1.Children(1).Visible='off';
    axis1.Children(2).Position=[350 50 10];
    axis1.Children(2).Color=[0.5 0.5 0.5];
else
    set(gca,'CameraPosition',[-1300 750 -25],'CameraTarget', [67 100 97],...
        'CameraUpVector', [0.11 0.04 -0.99],'CameraViewAngle', 2);
    axis1.Children(1).Visible='off';
    axis1.Children(2).Position=[-350 500 100];
    axis1.Children(2).Color=[0.5 0.5 0.5];
end

print(Fig, ['Subj_' SubjID '_' hs(ihs) 'cr_' Modality '_' Var '_A1-PT' Sufix '.tif'], '-dtiff')


%%
if ihs==1
    set(gca,'CameraPosition',[-1600 -1700 550],'CameraTarget',[140 115 89.5],...
        'CameraUpVector',[-0.13 -0.14 -0.98],'CameraViewAngle',1);
    axis1.Children(1).Visible='off';
    axis1.Children(2).Position=[0 -90 270];
    axis1.Children(2).Color=[0.5 0.5 0.5];
else
    set(gca,'CameraPosition',[1860 -1720 550],'CameraTarget', [76 57 106],...
        'CameraUpVector', [0.12 -0.12 -0.98],'CameraViewAngle', 1);
    axis1.Children(2).Position=[200 -200 100];
    axis1.Children(2).Color=[0.9 0.9 0.9];
end

print(Fig, ['Subj_' SubjID '_' hs(ihs) 'cr_' Modality '_' Var '_V1' Sufix '.tif'], '-dtiff')


%%
if ihs==1
    set(gca,'CameraPosition',[-451 -2250 975],'CameraTarget',[112 87 86.6],...
        'CameraUpVector',[-0.08 -0.33 -0.94],'CameraViewAngle',1);
    axis1.Children(1).Visible='off';
    axis1.Children(2).Position=[-200 -100 300];
    axis1.Children(2).Color=[0.65 0.65 0.65];
else
    set(gca,'CameraPosition',[720 -2190 1140],'CameraTarget', [73 55 103],...
        'CameraUpVector', [0.11 -0.3 -0.91],'CameraViewAngle', 1);
    axis1.Children(2).Position=[150 -150 500];
    axis1.Children(2).Color=[0.5 0.5 0.5];
end

print(Fig, ['Subj_' SubjID '_' hs(ihs) 'cr_' Modality '_' Var '_V2-3' Sufix '.tif'], '-dtiff')

end

