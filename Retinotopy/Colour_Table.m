%% Retinotopy
cycles = 2;

steps = 255 / 4 / cycles;   % steps from one colour to the next

% colour peaks
R = [1 0 0];
Y = [1 1 0];
G = [0 1 0];
B = [0 0 1];

% create the colour map
cmap = [];
for c = 1:cycles
    cmap = [cmap; linspace(R(1), Y(1), steps)', linspace(R(2), Y(2), steps)', linspace(R(3), Y(3), steps)'];
    cmap = [cmap; linspace(Y(1), G(1), steps)', linspace(Y(2), G(2), steps)', linspace(Y(3), G(3), steps)'];
    cmap = [cmap; linspace(G(1), B(1), steps)', linspace(G(2), B(2), steps)', linspace(G(3), B(3), steps)'];
    cmap = [cmap; linspace(B(1), R(1), steps)', linspace(B(2), R(2), steps)', linspace(B(3), R(3), steps)'];
end

img = repmat(cmap(:,1),1,100);
img(:,:,2) = repmat(cmap(:,2),1,100);
img(:,:,3) = repmat(cmap(:,3),1,100);

imagesc(img)

set(gca,'tickdir', 'out', ...
    'xtick', [], ...
    'xticklabel', [], ...
    'ytick', linspace(1,size(cmap(:,1),1),9), ...
    'yticklabel', linspace(0,360,9))

FileName = 'colortable-retinotopy-rgb.xml';

fid = fopen(FileName, 'w'); 

range = linspace(pi*-1, pi, size(cmap(:,1),1));

fprintf(fid, '<ColorMap name="retinotopy" space="RGB">\n');
fprintf(fid, '  <Point x="0" o="0" r="0.5" g="0.5" b="0.5"/>\n');
for i=1:size(cmap,1)
fprintf(fid, '  <Point x="%i" o="0.5" r="%f" g="%f" b="%f"/>\n', ...
    range(i), cmap(i,1), cmap(i,2), cmap(i,3));
end
fprintf(fid, '</ColorMap>');

fclose(fid);


%% MNI ROI
ROIs = {...
'Te10_L_MNI', 1, [1 0 0] ; ...
'Te11_L_MNI', 2, [1 0.2 0.2] ; ...
'Te12_L_MNI', 3, [1 0.4 0.4] ; ...
'TI1_L_MNI', 4, [1 0.6 0.6] ; ...

'hOC1_L_MNI', 11, [0 1 0] ; ...
'hOc2_L_MNI', 12, [0.2 1 0.2] ; ...
'hOc3d_L_MNI', 13, [0.4 1 0.4] ; ...
'hOc3v_L_MNI', 13.5, [0.4 1 0.4] ; ...
'hOc4d_L_MNI', 14, [0.6 1 0.6] ; ...
'hOc4v_L_MNI', 14.5, [0.6 1 0.6] ; ...
'hOc5_L_MNI', 15, [0.8 1 0.8] ; ...

'Te10_R_MNI', 101, [1 0 0] ; ...
'Te11_R_MNI', 102, [1 0.1 0.1] ; ...
'Te12_R_MNI', 103, [1 0.2 0.2] ; ...
'TI1_R_MNI', 104, [1 0.3 0.3] ; ...

'hOC1_R_MNI', 111, [0 1 0] ; ...
'hOc2_R_MNI', 112, [0.1 1 0.1] ; ...
'hOc3d_R_MNI', 113, [0.2 1 0.2] ; ...
'hOc3v_R_MNI', 113.5, [0.3 1 0.3] ; ...
'hOc4d_R_MNI', 114, [0.4 1 0.4] ; ...
'hOc4v_R_MNI', 114.5, [0.5 1 0.5] ; ...
'hOc5_R_MNI', 115, [0.6 1 0.6] ; ...
};

FileName = 'colortable-MNI_ROI-rgb.xml';

fid = fopen(FileName, 'w'); 

fprintf(fid, '<ColorMap name="MNI_ROI" space="RGB">\n');
fprintf(fid, '  <Point x="0" o="1" r="0.5" g="0.5" b="0.5"/>\n');
for i=1:size(ROIs,1)
fprintf(fid, '  <Point x="%i" o="1" r="%f" g="%f" b="%f"/>\n', ...
    ROIs{i,2}, ROIs{i,3}(1), ROIs{i,3}(2), ROIs{i,3}(3));
end
fprintf(fid, '</ColorMap>');

fclose(fid);


%% Heat map
A = colormap('hot')

img = repmat(A(:,1),1,100);
img(:,:,2) = repmat(A(:,2),1,100);
img(:,:,3) = repmat(A(:,3),1,100);

imagesc(img)

set(gca,'tickdir', 'out', ...
    'xtick', [], ...
    'xticklabel', [], ...
    'ytick', linspace(1,256,9), ...
    'yticklabel', linspace(0,360,9))

range = linspace(0, 4000, 256);

FileName = 'colortable-hot-rgb.xml';

fid = fopen(FileName, 'w'); 

fprintf(fid, '<ColorMap name="hot" space="RGB">\n');
% fprintf(fid, '  <Point x="0" o="1" r="0.5" g="0.5" b="0.5"/>\n');
for i=1:size(ROIs,1)
fprintf(fid, '  <Point x="%i" o="1" r="%f" g="%f" b="%f"/>\n', ...
    range(i), A(i,1),  A(i,2),  A(i,3));
end
fprintf(fid, '</ColorMap>');

fclose(fid);