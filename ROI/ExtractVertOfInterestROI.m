%%
clc; clear; close all

StartFolder=fullfile(pwd, '..','..');
addpath(genpath(fullfile(StartFolder, 'SubFun')))

SubjectList = [...
    '02';...
    '03';...
    '04';...
    '06';...
    '07';...
    '08';...
    '09';...
    '11';...
    '12';...
    '13';...
    '14';...
    '15';...
    '16'
    ];




for SubjInd = 1:size(SubjectList,1)
    
    fprintf('\n\n\n')
    
    % Subject's Identity
    SubjID = SubjectList(SubjInd,:);
    fprintf('Subject %s\n',SubjID)
    
    SubjectFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjID]);

    for hs=1:2
       
        if hs==1
            HsPrefix = 'l';
            fprintf('Left HS\n')
        else
             HsPrefix = 'r';
             fprintf('Right HS\n')
        end
        
        
        
        %% Get surface
        cd(fullfile(SubjectFolder, 'Structural', 'CBS'))

        vtk = dir(['T1_' SubjID '_thresh_clone_transform_strip_clone_transform_bound_mems_' ...
            HsPrefix 'cr_gm_avg.vtk']);
        [vertex,faces,~] = read_vtk(fullfile(SubjectFolder, 'Structural', 'CBS', vtk.name), 0, 1);

        NbVertex(hs) = size(vertex,2); %#ok<*SAGROW>
        
        
        
        %% A1
        cd(fullfile(SubjectFolder, 'ROI_MIPAV', 'A1'))

        vtk = dir(['Subj_' SubjID '_A1_' HsPrefix 'cr_norm_RG_UN.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('A1\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))
        
        ROI(1).name = 'A1';
        ROI(1).VertOfInt{hs} = find(mapping>0);
        clear mapping

        
        
        %% PT
        cd(fullfile(SubjectFolder, 'PMap'))

        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_PT_thres40.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('PT\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(2).name = 'PT';          
        ROI(2).VertOfInt{hs} = find(mapping==1);
        clear mapping



        %% V1 and V2-3
        cd(fullfile(SubjectFolder, 'PMap'))
       
        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_V1_V2-3_thres10.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('V1 - V2-3\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(3).name = 'V1';
        ROI(3).VertOfInt{hs} = find(mapping==1);

        ROI(4).name = 'V2-3';
        ROI(4).VertOfInt{hs} = find(mapping==2);
        clear mapping

        
        %% V1 and V2-3 Act/Deact
        cd(fullfile(SubjectFolder, 'ROI_MIPAV'))

        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_V1_act-deact.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('V1 act/deact\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(5).name = 'V1_act';
        ROI(5).VertOfInt{hs} = find(mapping==1);

        ROI(6).name = 'V1_deact';
        ROI(6).VertOfInt{hs} = find(mapping==-1);
        clear mapping
        
        
        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_V23_act-deact.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('V23 act/deact\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(7).name = 'V23_act';
        ROI(7).VertOfInt{hs} = find(mapping==1);

        ROI(8).name = 'V23_deact';
        ROI(8).VertOfInt{hs} = find(mapping==-1);
        clear mapping
        
        
        %% All ROIS Act/Deact for A and V
        
        Conditions_Names = {...
            'A Stim - Auditory Attention', ...
            'V Stim - Auditory Attention', ...
            'AV Stim - Auditory Attention', ...
            'A Stim - Visual Attention', ...
            'V Stim - Visual Attention', ...
            'AV Stim - Visual Attention'};
        
        ROI_idx = 9;
        for iCond=1:2
            for iROI=1:4
                vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_' ROI(iROI).name ...
                    '_' Conditions_Names{iCond}(1) '_act-deact.vtk']);
                [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
                
                fprintf('%s %s act/deact\n', ROI(iROI).name, Conditions_Names{iCond}(1))
                if numel(mapping)~=NbVertex(hs)
                    warning('The number of vertices does not match.')
                    numel(mapping)
                    NbVertex(hs)
                end
                disp(unique(mapping))
                
                ROI(ROI_idx).name = [ROI(iROI).name '_' Conditions_Names{iCond}(1) '_act'];
                ROI(ROI_idx).VertOfInt{hs} = find(mapping==1);
                ROI_idx = ROI_idx + 1;
                
                ROI(ROI_idx).name = [ROI(iROI).name '_' Conditions_Names{iCond}(1) '_deact'];
                ROI(ROI_idx).VertOfInt{hs} = find(mapping==-1);
                ROI_idx = ROI_idx + 1;
                clear mapping
            end
        end
        
        
        %% A1 and PT Act/Deact
        cd(fullfile(SubjectFolder, 'ROI_MIPAV'))

        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_A1_act-deact.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('A1 act/deact\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(25).name = 'A1_act';
        ROI(25).VertOfInt{hs} = find(mapping==1);

        ROI(26).name = 'A1_deact';
        ROI(26).VertOfInt{hs} = find(mapping==-1);
        clear mapping
        
        
        vtk = dir(['Subj_' SubjID '_' HsPrefix 'cr_PT_act-deact.vtk']);
        [~,~,mapping] = read_vtk(fullfile(pwd, vtk.name), 0, 1);
        
        fprintf('PT act/deact\n')
        if numel(mapping)~=NbVertex(hs)
            warning('The number of vertices does not match.')
            numel(mapping)
            NbVertex(hs)
        end
        disp(unique(mapping))

        ROI(27).name = 'PT_act';
        ROI(27).VertOfInt{hs} = find(mapping==1);

        ROI(28).name = 'PT_deact';
        ROI(28).VertOfInt{hs} = find(mapping==-1);
        clear mapping
        
        

        %%
        cd(fullfile(SubjectFolder, 'ROI_MIPAV'));
        mapping=zeros(1,size(vertex,2));
        mapping(ROI(1).VertOfInt{hs})=1;
        mapping(ROI(2).VertOfInt{hs})=2;
        mapping(ROI(3).VertOfInt{hs})=10;
        mapping(ROI(4).VertOfInt{hs})=15;
        write_vtk(['Subj_' SubjID '_' HsPrefix 'cr_AllROIs.vtk'], vertex, faces, mapping)
    
        
    end
    
    %%
    save(fullfile(SubjectFolder,'Transfer','ROI',['Subj_' SubjID '_ROI_VertOfInt.mat']), 'ROI', 'NbVertex')
    
    clear NbVertex ROI

end



cd(StartFolder)