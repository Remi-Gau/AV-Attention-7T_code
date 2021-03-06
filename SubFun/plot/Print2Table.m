function Print2Table(fid, ROIs, iROI, DATA, OneSideTTest)

ToPermute = DATA.ToPermute;

fprintf (fid, '%s,', ROIs{iROI});
for i=1:size(DATA.Betas,2)
    
    clear Perms
    
    for iPerm = 1:size(ToPermute,1)
        tmp2 = ToPermute(iPerm,:)';
        Perms(iPerm,:) = mean(DATA.Betas(:,i).*repmat(tmp2,1,size(DATA.Betas(:,i),2))); %#ok<AGROW,*SAGROW>
    end
    
    fprintf (fid, '%f,',nanmean(DATA.Betas(:,i)));
    fprintf (fid, '(,');
    fprintf (fid, '%f,',nanstd(DATA.Betas(:,i)));
    fprintf (fid, '),');
    
    if exist('OneSideTTest', 'var')
        if strcmp(OneSideTTest{i},'left')
            P(i) = sum(Perms<mean(DATA.Betas(:,i)))/numel(Perms); %#ok<*AGROW>
        elseif strcmp(OneSideTTest{i},'right')
            P(i)  = sum(Perms>mean(DATA.Betas(:,i)))/numel(Perms);
        elseif strcmp(OneSideTTest{i},'both')
            P(i)  = sum(abs(Perms)>abs(mean(DATA.Betas(:,i)))) / numel(Perms) ;
        end
        %     [~,P(i),~,STATS] = ttest(DATA.Betas(:,i), 0, 0.05, OneSideTTest{i});
        
    else
        P(i)  = sum(abs(Perms)>abs(mean(DATA.Betas(:,i)))) / numel(Perms) ;
        %     [~,P(i),~,STATS] = ttest(DATA.Betas(:,i), 0, 0.05);
    end
    if (exist('OneSideTTest', 'var') && ~strcmp(OneSideTTest{i},'both'))
        fprintf (fid, '1 sided ');
    else
        fprintf (fid, '2 sided ');
    end
    
    %     fprintf (fid, '%f,',STATS.tstat);
    fprintf (fid, ',');
    if P(i)<0.001
        fprintf (fid, '<.001,');
    else
    fprintf (fid, '%f,',P(i));
    end
    
    % Print effect size
    CI = bootci(1000,{@(x) Unbiased_ES(x), DATA.Betas(:,i)},'alpha',0.05,'type','bca');
    fprintf (fid, '[,%f,-,%f,]',CI(1),CI(2));
    fprintf (fid, ',');
  
end

fprintf (fid, '\n');


end

function du = Unbiased_ES(grp_data)
% from DOI 10.1177/0013164404264850
d = mean(grp_data)/std(grp_data);
nu = length(grp_data)-1;
G = gamma(nu/2)/(sqrt(nu/2)*gamma((nu-1)/2));
du = d*G;
end