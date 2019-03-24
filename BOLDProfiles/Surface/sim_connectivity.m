clear all; close all;

num_sim = 1000;

mean_act = 2;
mean_deact = -0.5;
std_act = 2;
std_deact = 1;

for kk = 1 : num_sim

    n = 12; % number of blocks for A and V stimulation
    num_vox_A = 1000;  % number of voxels in A1
    num_vox_V = 3000;  % number of voxels in A1
    
    %------------------------------------------------------------
    % generate 'true' signal for A and V stim conditions in A1
    mu_A_A1 = repmat(mean_act,n,num_vox_A); 
    mu_V_A1 = repmat(mean_deact,n,num_vox_A);

    % Add noise
    x_A_A1 = mu_A_A1 + randn(n,num_vox_A)*std_act;
    x_V_A1 = mu_V_A1 + randn(n,num_vox_A)*std_deact;
    
    % Response over both stim conditions
    x_A1 = [x_A_A1; x_V_A1];
    
    %------------------------------------------------------------
    % generate 'true' signal for A and V stim conditions in V1
    mu_A_V1 = repmat(mean_deact,n,num_vox_A); 
    mu_V_V1 = repmat(mean_act,n,num_vox_A);

    % Add noise
    x_A_V1 = mu_A_A1 + randn(n,num_vox_A)*std_deact;
    x_V_V1 = mu_V_A1 + randn(n,num_vox_A)*std_act;
    
    % Response over both stim conditions
    x_V1 = [x_A_V1; x_V_V1];
    
    
    
    
    %------------------------------------------------------------
    % Regress V activity against A activations separately for voxels Astim
    % >0 and Astim<=0
    for my_pos = 1 : size(x_A1_pos,2)
        beta_pos(my_pos) = regress(y_V1,x_A1_pos(:,my_pos));
    end

    for my_neg = 1 : size(x_A1_neg,2)
        beta_neg(my_neg) = regress(y_V1,x_A1_neg(:,my_neg));
    end

    mean_pos(kk) = mean(beta_pos);
    mean_neg(kk) =mean(beta_neg);
    
end

figure
plot(mean_pos,'b:'); hold on; plot(mean_neg,'r:')

