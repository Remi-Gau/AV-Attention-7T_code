function PrintTableCorrCoeff(SavedTxt, ROI, ToPlot, Effect, Values2Plot,ToPermute)

if isempty(ToPermute)
    [~,P,~,STATS] = ttest(Values2Plot);
else
    Perms = ToPermute.*repmat(Values2Plot,[size(ToPermute,1),1]);
    Perms = mean(Perms,2);
    P = sum( abs(Perms) > abs( mean(Values2Plot) ) ) / numel(Perms);
end

Legends = {'mean', '(','STD',')', 't value','p value', 'effect size'};

fid = fopen (SavedTxt, 'w');
fprintf (fid, 'Correlation coefficient between profiles\n');
fprintf (fid, '%s,%s,%s\n', ROI, ToPlot, Effect);
for i=1:length(Legends)
    fprintf (fid, '%s,', Legends{i});
end
fprintf (fid, '\n');

fprintf (fid, '%f,',nanmean(Values2Plot));

fprintf (fid, '(,');
fprintf (fid, '%f,',nanstd(Values2Plot));
fprintf (fid, '),');

fprintf (fid, '2 sided ');
if isempty(ToPermute)
    fprintf (fid, '%f',STATS.tstat);
end
fprintf (fid, ',');

if P<0.001
    fprintf (fid, '<.001,');
else
    fprintf (fid, '%f,',P);
end

%     fprintf (fid, '%f,',abs(nanmean(Values2Plot)/nanstd(Values2Plot)));

    % Print effect size
    CI = bootci(1000,{@(x) Unbiased_ES(x), Values2Plot},'alpha',0.05,'type','bca');
    fprintf (fid, '[,%f,-,%f,]',CI(1),CI(2));
    fprintf (fid, ',');
 
    fclose (fid);

end


function du = Unbiased_ES(grp_data)
% from DOI 10.1177/0013164404264850
d = mean(grp_data)/std(grp_data);
nu = length(grp_data)-1;
G = gamma(nu/2)/(sqrt(nu/2)*gamma((nu-1)/2));
du = d*G;
end