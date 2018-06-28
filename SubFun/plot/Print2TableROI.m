function Print2TableROI(fid, ROIs, iROI, DATA, OneSideTTest)

if DATA.MVPA
    DATA.Data = DATA.Data-0.5;
end

ToPermute = DATA.ToPermute;
for iPerm = 1:size(ToPermute,1)
    tmp2 = ToPermute(iPerm,:)';
    Perms(iPerm,:) = mean(DATA.Data.*repmat(tmp2,1,size(DATA.Data,2))); %#ok<*SAGROW>
end

fprintf (fid, '%s,', ROIs{iROI});
fprintf (fid, '%f,',nanmean(DATA.Data));
fprintf (fid, '(,');
fprintf (fid, '%f,',nanstd(DATA.Data));
fprintf (fid, '),');

% if exist('OneSideTTest', 'var')
%     [~,P] = ttest(DATA.Data, 0, 'alpha', 0.05, 'tail', OneSideTTest);
% else
%     [~,P] = ttest(DATA.Data, 0, 'alpha', 0.05);
% end
if isfield(DATA, 'OneSideTTest')
    if strcmp(DATA.OneSideTTest,'left')
        P = sum(Perms<mean(DATA.Data))/numel(Perms);
    elseif strcmp(DATA.OneSideTTest,'right')
        P = sum(Perms>mean(DATA.Data))/numel(Perms);
    elseif strcmp(DATA.OneSideTTest,'both')
        P  = sum(abs(Perms)>abs(mean(DATA.Data))) / numel(Perms) ;
    end
else
    P  = sum(abs(Perms)>abs(mean(DATA.Data))) / numel(Perms) ;
end

    if (isfield(DATA, 'OneSideTTest') && ~strcmp(DATA.OneSideTTest,'both'))
        fprintf (fid, '1 sided ');
    else
        fprintf (fid, '2 sided ');
    end
    
    %     fprintf (fid, '%f,',STATS.tstat);
    fprintf (fid, ',');
    if P<0.001
        fprintf (fid, '<.001,');
    else
    fprintf (fid, '%f,',P);
    end
    
fprintf (fid, '%f,\n',abs(nanmean(DATA.Data)/nanstd(DATA.Data)));


end

