function [P,H] = SignPermTest(tmp)
sets = {};
for iSub=1:11
    sets{iSub} = [-1 1]; %#ok<*AGROW>
end
[a, b, c, d, e, f, g, h, i, j, k] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:)];


for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:);
    tmp2 = repmat(tmp2',1,size(tmp,2));
    Perms(iPerm,:) = mean(tmp.*tmp2);  %#ok<*AGROW>
end
P = sum( ...
    abs( Perms ) > ...
    repmat( abs(mean(tmp)), size(Perms,1),1)  ) ...
    / size(Perms,1);

if P<.05
    H = 1;
else
    H = 0;
end
end