function PlotRectangle(NbLayers,Fontsize,Switch,LabelDepth)

if nargin<4 || isempty(LabelDepth)
    LabelDepth=1;
end

switch NbLayers
    case 8
        COLOR_Layer= flipud([
            255,247,236;
            254,232,200;
            253,212,158;
            253,187,132;
            252,141,89;
            239,101,72;
            215,48,31;
            153,0,0]);

    case 6
        COLOR_Layer= flipud([
            254,229,217;
            252,187,161;
            252,146,114;
            251,106,74;
            222,45,38;
            165,15,21]);
%         COLOR_Layer = repmat(linspace(.7,.3,6)',1,3)*255;
    case 4
        COLOR_Layer= [
            254,229,217;
            252,174,145;
            251,106,74;
            203,24,29];
end

COLOR_Layer = COLOR_Layer/255;

ax = gca;
axPos = ax.Position;
axPos(2) = axPos(2)-.03;
axPos(4) = .02;
axes('Position',axPos);
% axis('off')
TEXT = round(linspace(0,100,NbLayers+2));
TEXT([1 end]) = [];
if Switch
    TEXT=fliplr(TEXT);    
end
if NbLayers == 8
    TEXT = round(linspace(0,100,NbLayers));
    TEXT=fliplr(TEXT);
end

RecPos = linspace(0,0.9,NbLayers+1);

for i=1:size(COLOR_Layer,1)
    rectangle('Position', [RecPos(i) 0 diff(RecPos(1:2)) 1], 'facecolor', COLOR_Layer(i,:), 'edgecolor', 'w');
    if LabelDepth
        t = text(RecPos(i)+diff(RecPos(1:2))/2-.023,0.5,num2str(TEXT(i)));
        set(t,'fontsize',Fontsize-4);
    end
end
axis([0 0.9 0 1])

set(gca,'color', 'none', 'tickdir', 'out', 'xtick', [0 0.45 .9],'xticklabel',  {'WM      ' 'GM' '      CSF'}, ...
    'ytick', [],'yticklabel', [], ...
    'ticklength', [0.00001 0], 'fontsize', 7)


end

