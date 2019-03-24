clc; clear;

V1_2_A1 = 1;
ROI_Suffix = 'V1_2_A1';

NbLayers = 7;
NbLayers = NbLayers+1;

load(fullfile(pwd,strcat('Connectivity_all_', ROI_Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')))

% Mean for each condition and block for all subjects (72 rows X 11 columns)
Grp_Cst_V1 = Grp_Cst_Seed;
Grp_Cst_A1 = Grp_Cst_Target;

% 1rst eingenvariate for each condition and block for all subjects (72 rows X 11 columns)
Grp_U_V1 = Grp_U_Seed;
Grp_U_A1 = Grp_U_Target;

% Vector defining which condition corresponds to which row
CdtVec = repmat(1:6,12,1);
CdtVec = CdtVec(:);

A = repmat([0 1;0 1],1,6);
A = logical(A(:));

% loop over subjects
for iSubj = 1:size(Grp_Cst_A1,2)
    
    close all
    
    % extract data from A and V stim under A att
    Row2Select = ismember(CdtVec,1:2)';
    X_V1_A_att = Grp_Cst_V1(Row2Select,iSubj);
    Y_A1_A_att = Grp_Cst_A1(Row2Select,iSubj);
    
    % extract data from A and V stim under V att
    Row2Select = ismember(CdtVec,4:5)';
    X_V1_V_att = Grp_Cst_V1(Row2Select,iSubj);
    Y_A1_V_att = Grp_Cst_A1(Row2Select,iSubj);
    
    Connectivity_vectors{iSubj} = [...
        X_V1_A_att Y_A1_A_att X_V1_V_att Y_A1_V_att];

    X_V1_A_att = Connectivity_vectors{iSubj}(:,1);
    Y_A1_A_att = Connectivity_vectors{iSubj}(:,2);
    X_V1_V_att = Connectivity_vectors{iSubj}(:,3);
    Y_A1_V_att = Connectivity_vectors{iSubj}(:,4);
    
    
    figure(1)
    
    subplot (121)
    
    scatter(X_V1_A_att, Y_A1_A_att,'r'); hold on
    scatter(X_V1_A_att(A), Y_A1_A_att(A),'.r'); hold on
    [B1,BINT1,R1] = regress(X_V1_A_att,  [Y_A1_A_att ones(24,1)])
    x = [-2:0.01:4]'
    plot(x,[x*B1(1)+ B1(2)*ones(length(x),1)],'r')
    
    scatter(X_V1_V_att, Y_A1_V_att,'b')
    scatter(X_V1_V_att(A), Y_A1_V_att(A),'.b')
    [B2,BINT2,R2] = regress(X_V1_V_att,  [Y_A1_V_att ones(24,1)])
    x = [-2:0.01:4]'
    plot(x,[x*B2(1)+ B2(2)*ones(length(x),1)],'b')
    title('A1 to V1')
    

    subplot (122)
    
    scatter(Y_A1_A_att,X_V1_A_att, 'r'); hold on
    scatter(Y_A1_A_att(A),X_V1_A_att(A), '.r'); hold on
    [B1,BINT1,R1] = regress(Y_A1_A_att ,  [X_V1_A_att ones(24,1)])
    x = [-2:0.01:4]'
    plot(x,[x*B1(1)+ B1(2)*ones(length(x),1)],'r')
    
    scatter(Y_A1_V_att,X_V1_V_att, 'b')
    scatter(Y_A1_V_att(A),X_V1_V_att(A), '.b')
    [B2,BINT2,R2] = regress(Y_A1_V_att ,  [X_V1_V_att ones(24,1)])
    x = [-2:0.01:4]'
    plot(x,[x*B2(1)+ B2(2)*ones(length(x),1)],'b')
    title('V1 to A1')
    
    
    pause(1)
end


save('Connectivity_vectors.mat', 'Connectivity_vectors')