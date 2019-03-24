clear
clc
close all

% A stim
% A1		1.580		(	0.440	)
% V1		-0.869		(	0.441	)
% V stim						
% A1		-0.096		(	0.154	)
% V1		1.467		(	0.597	)

% set group mean and STD
% mean_act_A1 = 1.58;
% mean_deact_A1 = -0.1;
% 
% mean_act_V1 = 1.47;
% mean_deact_V1 = -0.87;
% 
% grp_std_act_A1 =.44;
% grp_std_deact_A1 =.15;
% 
% grp_std_act_V1 =.59;
% grp_std_deact_V1 =.44;

mean_act_A1 = 1.58;
mean_deact_A1 = -0.1;

mean_act_V1 = 1.47;
mean_deact_V1 = -0.87;

grp_std_act_A1 =.44;
grp_std_deact_A1 =.15;

grp_std_act_V1 =.59;
grp_std_deact_V1 =.44;

% set within subject parameters
std_act_A1 = 4;
std_deact_A1 = 2;

std_act_V1 = 8;
std_deact_V1 = 4;

% set attention effect
att_diff_A1 = [-0.18 -0.06];
att_diff_V1 = [-0.01 0.12];
% att_diff_A1 = [0 0];
% att_diff_V1 = [0 0];

num_subj = 3;

n = 12; % number of blocks for A and V stimulation

num_vox_A = 500;  % number of voxels in A1
num_vox_V = num_vox_A*3;  % number of voxels in V1

for kk = 1 : num_subj
    
    fprintf('Subject %i/%i\n',kk,num_subj)
    
    if kk==1
        figure(1)
    end
    
    % set subject specific noise
    Subj_spe_effect_act_A1 = randn(1)*grp_std_act_A1;
    Subj_spe_effect_deact_A1 = randn(1)*grp_std_deact_A1;
    Subj_spe_effect_act_V1 = randn(1)*grp_std_act_V1;
    Subj_spe_effect_deact_V1 = randn(1)*grp_std_deact_V1;
    
    for iAtt = 1:2
        
        
        %% A1
        % generate 'true' signal for A and V stim conditions in A1
        mu_A_A1 = repmat(mean_act_A1+Subj_spe_effect_act_A1,n,num_vox_A);
        mu_V_A1 = repmat(mean_deact_A1+Subj_spe_effect_deact_A1,n,num_vox_A);
        
        % Add noise
        x_A_A1 = mu_A_A1 + randn(n,num_vox_A)*std_act_A1;
        x_V_A1 = mu_V_A1 + randn(n,num_vox_A)*std_deact_A1;
        
        % Response over both stim conditions
        x_A1 = [x_A_A1; x_V_A1];
        
        
        %% V1
        % generate 'true' signal for A and V stim conditions in V1
        mu_A_V1 = repmat(mean_deact_V1+Subj_spe_effect_deact_V1,n,num_vox_V);
        mu_V_V1 = repmat(mean_act_V1+Subj_spe_effect_act_V1,n,num_vox_V);
        
        % Add noise
        x_A_V1 = mu_A_V1 + randn(n,num_vox_V)*std_deact_V1;
        x_V_V1 = mu_V_V1 + randn(n,num_vox_V)*std_act_V1;
        
        % Response over both stim conditions
        x_V1 = [x_A_V1; x_V_V1];
        
        
        %% Attention modulation
        if iAtt==2
            x_V1 = x_V1+[ones(n,1)*att_diff_V1(1);ones(n,1)*att_diff_V1(2)];
            x_A1 = x_A1+[ones(n,1)*att_diff_A1(1);ones(n,1)*att_diff_A1(2)];
        end
        
        
        %% Plot
        if kk==1
            if iAtt==1
                Subplot = 1:4;
            else
                Subplot = 5:8;
            end
            
            for i = 1:4
                switch i
                    case 1
                        Data2Plot = x_V1(1:n,:);
                    case 2
                        Data2Plot = x_A1(1:n,:);
                    case 3
                        Data2Plot = x_V1((n+1):end,:);
                    case 4
                        Data2Plot = x_A1((n+1):end,:);
                end
                
                subplot(4,2,Subplot(i))
                
                hold on
                hist(mean(Data2Plot),100)
                
                axis tight
                ax = axis;
                plot([0 0], [0 ax(4)], 'k', 'linewidth',2)
                plot([mean(mean(Data2Plot)) mean(mean(Data2Plot))], [0 ax(4)], '-r',...
                    'linewidth',2)
                axis([-8 8 0 ax(4)])
            end
        end
        
        
        %% Do regression after averaging
        X = mean(x_V1,2);
        Y = mean(x_A1,2);
        B(:,kk,1,iAtt) = glmfit(X,Y,'normal'); % Regression V1 to A1
        
        
        %% Do regression before averaging
        % V1 to A1
        for iVert = 1:size(x_V1,2) % regress for each vertex
            tmp(:,iVert) = glmfit(...
                x_V1(:,iVert),...
                Y,'normal');
        end
        B(:,kk,2,iAtt)=mean(tmp,2); % average
        
    end
    
end

%% Add labels
subplot(4,2,1)
title('V1')
ylabel('Stim A_{att A}')
subplot(4,2,2)
title('A1')
subplot(4,2,3)
ylabel('Stim V_{att A}')
subplot(4,2,5)
ylabel('Stim A_{att V}')
subplot(4,2,7)
ylabel('Stim V_{att V}')


%%
clc
fprintf('Connectivity: AVERAGE then REGRESS\n')
fprintf('Connectivity V1 to A1 for A att: %f +/- %f\n', ...
    mean(B(2,:,1,1)), std(B(2,:,1,1)))
fprintf('Connectivity V1 to A1 for V att: %f +/- %f\n', ...
    mean(B(2,:,1,2)), std(B(2,:,1,2)))
[~,P] = ttest(B(2,:,1,1)', B(2,:,1,2)')

fprintf('\n\n\n')
fprintf('Connectivity: REGRESS then AVERAGE\n')
fprintf('Connectivity V1 to A1 for A att: %f +/- %f\n', ...
    mean(B(2,:,2,1)), std(B(2,:,2,1)))
fprintf('Connectivity V1 to A1 for V att: %f +/- %f\n', ...
    mean(B(2,:,2,2)), std(B(2,:,2,2)))
[~,P] = ttest(B(2,:,2,1)', B(2,:,2,2)')


%% Bootstrap and p curves
NbSubj = 11;
for iBt=1:20000
    idx = randperm(num_subj);
    idx = idx(1:NbSubj);
%     P1(iBt) = SignPermTest(B(2,idx,1,1)'-B(2,idx,1,2)');
%     P1(iBt) = SignPermTest(B(2,idx,2,1)'-B(2,idx,2,2)');
    [H1(iBt),P1(iBt)] = ttest(B(2,idx,1,1)', B(2,idx,1,2)');
    [H2(iBt),P2(iBt)] = ttest(B(2,idx,2,1)', B(2,idx,2,2)');
end

%% Plot
close all

figure('name', 'p curves connectivity', 'Position', [50, 50, 1000, 700])
subplot(211)
hist(P1,1000)
title('Connectivity: AVERAGE then REGRESS')
xlabel('p value')
subplot(212)
hist(P2,1000)
title('Connectivity: REGRESS then AVERAGE')
xlabel('p value')


subplot(211)
ax=axis;
axis([0 1 0 ax(4)])
% ax=axis;
text(.65,ax(4)*.95,sprintf('act A1 = %.3f +/- %.3f',mean_act_A1,grp_std_act_A1))
text(.65,ax(4)*.90,sprintf('deact A1 = %.3f +/- %.3f',mean_deact_A1,grp_std_deact_A1))

text(.65,ax(4)*.85,sprintf('act V1 = %.3f +/- %.3f',mean_act_V1,grp_std_act_V1))
text(.65,ax(4)*.80,sprintf('deact V1 = %.3f +/- %.3f',mean_deact_V1,grp_std_deact_V1))


text(.65,ax(4)*.70,sprintf('w/in subj std act/deact A1 = %.3f / %.3f',std_act_A1,std_deact_A1))
text(.65,ax(4)*.65,sprintf('w/in subj std act/deact V1 = %.3f / %.3f',std_act_V1,std_deact_V1))

text(.60,ax(4)*.55,sprintf('attention mod A1: A stim / V stim= %.3f / %.3f',att_diff_A1(1),att_diff_A1(2)))
text(.60,ax(4)*.50,sprintf('attention mod V1: A stim / V stim= %.3f / %.3f',att_diff_V1(1),att_diff_V1(2)))


text(.1,ax(4)*.90,sprintf('Connectivity V1 to A1 for A att: %.3f +/- %.3f\n', ...
    mean(B(2,:,1,1)), std(B(2,:,1,1))))
text(.1,ax(4)*.85,sprintf('Connectivity V1 to A1 for V att: %.3f +/- %.3f\n', ...
    mean(B(2,:,1,2)), std(B(2,:,1,2))))

text(.1,ax(4)*.65,...
    sprintf('p( p_{regress then average}<.05 ) = %.3f',...
    sum(H1)/numel(H1)))


subplot(212)
ax=axis;
axis([0 1 0 ax(4)])

text(.1,ax(4)*.65,...
    sprintf('p( p_{average then regress}<.05 ) = %.3f',...
    sum(H2)/numel(H1)))

text(.1,ax(4)*.55,...
    sprintf('p( p_{regress then average}<.05 | p_{average then regress}>.05 ) = %.3f',...
    sum(H2>H1)/numel(H1)))

text(.1,ax(4)*.90,sprintf('Connectivity V1 to A1 for A att: %.3f +/- %.3f\n', ...
    mean(B(2,:,2,1)), std(B(2,:,2,1))))
text(.1,ax(4)*.85,sprintf('Connectivity V1 to A1 for V att: %.3f +/- %.3f\n', ...
    mean(B(2,:,2,2)), std(B(2,:,2,2))))