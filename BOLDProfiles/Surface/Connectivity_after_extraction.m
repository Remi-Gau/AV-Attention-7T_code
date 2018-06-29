clc; clear;

StartFolder=fullfile(pwd,'..','..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))
Get_dependencies('D:\Dropbox\')

COLOR_Subject = ColorSubject();

V1_2_A1 = 1;
ROI_Suffix = 'V1_2_A1';

NbLayers = 7;
NbLayers = NbLayers+1;

FigureFolder='D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Figures\Profiles\Surfaces\Baseline';
Results_Folder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\Results\Profiles\Surfaces';

load(fullfile(Results_Folder,strcat('Connectivity_all_', ROI_Suffix, '_Surf_', num2str(NbLayers), '_layers.mat')))


Grp_Cst_V1 = Grp_Cst_Seed;
Grp_Cst_A1 = Grp_Cst_Target;

Grp_U_V1 = Grp_U_Seed;
Grp_U_A1 = Grp_U_Target;

CdtVec = repmat(1:6,12,1);
CdtVec = CdtVec(:);

for iSubj = 1:size(Grp_Cst_A1,2)
    for i=1:2
        if i==1
            Row2Select = ismember(CdtVec,1:2)'; % for A and V stim under A att
        else
            Row2Select = ismember(CdtVec,4:5)'; % for A and V stim under V att
        end
        
%         X = Grp_Cst_V1(Row2Select,iSubj);
%         Y = Grp_Cst_A1(Row2Select,iSubj);
        
        X = Grp_U_Seed(Row2Select,iSubj);
        Y = Grp_U_A1(Row2Select,iSubj);
        
        Connectivity_A1_to_V1(iSubj,i,:) = regress(X,[Y ones(size(Y))]);
        Connectivity_V1_to_A1(iSubj,i,:) = regress(Y,[X ones(size(X))]);
        
        % No constant
%         Connectivity_A1_to_V1_no_cst(iSubj,i,:) = regress(X,Y);
%         Connectivity_V1_to_A1_no_cst(iSubj,i,:) = regress(Y,X);
        
    end
end


%%
clc

fprintf('Connectivity V1 to A1 for A att: %f ± %f\n', mean(Connectivity_V1_to_A1(:,1,1)) , nansem(Connectivity_V1_to_A1(:,1,1)))
P = SignPermTest(Connectivity_V1_to_A1(:,1,1))
fprintf('Connectivity V1 to A1 for V att: %f ± %f\n', mean(Connectivity_V1_to_A1(:,2,1)) , nansem(Connectivity_V1_to_A1(:,2,1)))
P = SignPermTest(Connectivity_V1_to_A1(:,2,1))
P = SignPermTest(Connectivity_V1_to_A1(:,1,1)-Connectivity_V1_to_A1(:,2,1))


fprintf('Connectivity A1 to V1 for A att: %f ± %f\n', mean(Connectivity_A1_to_V1(:,1,1)) , nansem(Connectivity_A1_to_V1(:,1,1)))
P = SignPermTest(Connectivity_A1_to_V1(:,1,1))
fprintf('Connectivity A1 to V1 for V att: %f ± %f\n', mean(Connectivity_A1_to_V1(:,2,1)) , nansem(Connectivity_A1_to_V1(:,2,1)))
P = SignPermTest(Connectivity_A1_to_V1(:,2,1))
P = SignPermTest(Connectivity_A1_to_V1(:,1,1)-Connectivity_A1_to_V1(:,2,1))


