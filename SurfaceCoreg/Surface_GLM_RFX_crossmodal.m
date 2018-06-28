%%
clear
clc

HS = 'L';
hs = 'l';
inf_hs = 'r';

StartFolder = fullfile(pwd, '..');
cd(StartFolder)
addpath(genpath(fullfile(StartFolder, 'SubFun')));
Get_dependencies('D:\Dropbox')

SubjToInclude = true(13,1);
SubjToInclude([4 11],1) = false;

NbLayers = 6;

Conditions_Names = {...
    'A-Stim_Auditory-Attention', ...
    'V-Stim_Auditory-Attention', ...
    'AV-Stim_Auditory-Attention', ...
    'A-Stim_Visual-Attention', ...
    'V-Stim_Visual-Attention', ...
    'AV-Stim_Visual-Attention'};

DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
% DesMat = [ones(NbLayers,1) DesMat'];
DesMat = spm_orth(DesMat);

Suffix = {'cst', 'lin', 'quad'};

DataFolder = 'D:\Dropbox\PhD\Experiments\AV_Integration_7T\archives\SurfCoreg';
cd(DataFolder)

for Smooth = 1
    
    
    InfSurfFile = fullfile(DataFolder, ['Surface_ls_' inf_hs 'h_inf.vtk']);
    [InfVertex,InfFace,InfMapping] = read_vtk(InfSurfFile, 0, 1);
    
    if Smooth
        suffix = '_smoothdata';
    else
        suffix = '';
    end
    
    %%
    clear AllLayers AllLayersNoSmooth
    
    for iCond = 1:6
        
        fprintf('\n Reading data\n')
        
        for iLayer = 1:NbLayers
            
            VTK_file = dir(['GrpSurf_' Conditions_Names{iCond} '_' num2str(iLayer) '-Layer_' hs 'h' suffix '.vtk']);
            
            disp(VTK_file.name)
            
            [Vertex,Face,Mapping] = read_vtk(VTK_file.name, 12, 1);
            
            AllLayers(:,:,iLayer,iCond) = Mapping;
            
            if Smooth
                VTK_file = dir(['GrpSurf_' Conditions_Names{iCond} '_' num2str(iLayer) '-Layer_' hs 'h.vtk']);
                [~,~,Mapping] = read_vtk(VTK_file.name, 12, 1);
                AllLayersNoSmooth(:,:,iLayer,iCond) = Mapping;
            end
            
        end
        
    end
    
    
    %% Cross modal
    fprintf('\n Cross modal')
    mkdir(fullfile(DataFolder, [HS 'H'], 'CrossModal'))
    CrossModal = {[3 6 ; 1 4] , [3 6 ; 2 5], [3;1], [6;4], [3;2], [6;5],};
    CrossModalName = {'AV-A', 'AV-V', 'AV-A_AttA', 'AV-A_AttV', 'AV-V_AttA', 'AV-V_AttV'};
    
    for iCrossModal = 1:numel(CrossModal)
        
        CrossModalLayers = ...
            sum(AllLayers(SubjToInclude,:,:,CrossModal{iCrossModal}(1,:)),4) ...
            - sum(AllLayers(SubjToInclude,:,:,CrossModal{iCrossModal}(2,:)),4);
        
        if Smooth
            CrossModalLayersNoSmooth = ...
                sum(AllLayersNoSmooth(SubjToInclude,:,:,CrossModal{iCrossModal}(1,:)),4) ...
                - sum(AllLayersNoSmooth(SubjToInclude,:,:,CrossModal{iCrossModal}(2,:)),4);
            IsZero = CrossModalLayersNoSmooth==0;
        else
            IsZero = CrossModalLayers==0;
        end
        IsZero = any(IsZero,3);
        tabulate(sum(IsZero))
        
        Cst = zeros(1, size(Mapping,2));
        Lin = zeros(1, size(Mapping,2));
        Quad = zeros(1, size(Mapping,2));
        
        for NbSub2Excl = 0:(sum(SubjToInclude)-4)
            
            fprintf('\nGLM on vertices with %i subjects', sum(SubjToInclude)-NbSub2Excl)
            
            VertOfInt = find(sum(IsZero)==NbSub2Excl);
            
            Y = [];
            for iVert = 1:numel(VertOfInt)
                Subj2Exclu = find(~IsZero(:,VertOfInt(iVert)));
                Y(:,iVert,:) = CrossModalLayers(Subj2Exclu,VertOfInt(iVert),:);
            end
            
            Y = shiftdim(Y,2);
            Y = reshape(Y, [size(Y,1)*size(Y,2),size(Y,3)] );
            
            X = [];
            for iSubj=1:(sum(SubjToInclude)-NbSub2Excl)
                X((1:6)+6*(iSubj-1),(1:size(DesMat,2))+size(DesMat,2)*(iSubj-1)) = DesMat; %#ok<SAGROW>
            end
            
            B = pinv(X)*Y;
            
            Cst_tmp = mean(B(1:3:size(X,2),:));
            Cst(VertOfInt) = Cst_tmp;
            
            Lin_tmp = mean(B(2:3:size(X,2),:));
            Lin(VertOfInt) = Lin_tmp;
            
            Quad_tmp = mean(B(3:3:size(X,2),:));
            Quad(VertOfInt) = Quad_tmp;
 
        end
        
        fprintf('\n');
        
                    write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
                        [CrossModalName{iCrossModal} '_' hs 'h_lin_inf' suffix '.vtk']), InfVertex, InfFace, Lin')
                    write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
                        [CrossModalName{iCrossModal} '_' hs 'h_cst_inf' suffix '.vtk']), InfVertex, InfFace, Cst')
        %             write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
        %                 [CrossModalName{iCrossModal} '_' hs 'h_quad_inf' suffix '.vtk']), InfVertex, InfFace, Quad')
        %
                    write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
                        [CrossModalName{iCrossModal} '_' hs 'h_lin' suffix '.vtk']), Vertex, Face, Lin')
                    write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
                        [CrossModalName{iCrossModal} '_' hs 'h_cst' suffix '.vtk']), Vertex, Face, Cst')
        %             write_vtk(fullfile(DataFolder, [HS 'H'], 'CrossModal', ...
        %                 [CrossModalName{iCrossModal} '_' hs 'h_quad' suffix '.vtk']), Vertex, Face, Quad')
        
    end
    
end