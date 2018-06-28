function ReadFeatVTK(AnalysisFolder, NbLayers, ROIs, NbVertices, BetaList, Betas, FilesCell, NbWorkers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

MatlabVer = version('-release');
if str2double(MatlabVer(1:4))>2013
    if isempty(gcp)
        KillGcpOnExit = 1;
        parpool(NbWorkers);
    else
        KillGcpOnExit = 0;
    end
else
    if matlabpool('size') == 0
        KillGcpOnExit = 1;
        matlabpool(NbWorkers)
    elseif matlabpool('size') ~= NbWorkers
        matlabpool close
        matlabpool(NbWorkers)
        KillGcpOnExit = 0;
    else
        KillGcpOnExit = 0;
    end
end

% Initialize variables
FeaturesAll = cell(size(ROIs,1),2,numel(NbLayers));
for iLayer = 1:numel(NbLayers)
    for hs = 1:2
        for iROI=1:size(ROIs,1)
            FeaturesAll{iROI,hs,iLayer} = nan(size(BetaList,1), numel(ROIs{iROI,hs+2})*NbLayers(iLayer));
        end
    end
end

% Gets the value for each ROI from the mapping of each beta image
for iLayer = 1:numel(NbLayers) % the mapping depends on the number of layers
    
    fprintf('  For %i layers\n', NbLayers(iLayer))
    
    % Format for reading the vertices from the VTK file
    Spec = repmat('%f ', 1, NbLayers(iLayer));
    
    % For the 2 hemispheres
    for hs = 1:2
        if hs==1
            fprintf('   Left hemipshere\n')
        else
            fprintf('   Right hemipshere\n')
        end
        
        AllMapping = nan(NbVertices(hs),NbLayers(iLayer),size(Betas,1));
        
        fprintf(1,'   [%s]\n   [ ',repmat('.',1,size(Betas,1)));
        
        parfor iBeta = 1:size(Betas,1)
            
            A = fileread(FilesCell{iBeta,hs,iLayer}); % reads file quickly
            B = A(strfind(A, 'TABLE default')+14:end); %clear A; % extracts lines that correspond to the mapping
            
            C = textscan(B, Spec, 'returnOnError', 0); %clear B; % extracts values from those lines
            mapping = cell2mat(C); %clear C
            
            if size(mapping,1)~=(NbVertices(hs)) %#ok<PFBNS>
                error('A VTK file has wrong number of vertices:\n%s', FilesCell{iBeta,hs,iLayer})
            end
            
            AllMapping(:,:,iBeta) = mapping;

            fprintf(1,'\b.\n');
            
        end      
        fprintf(1,'\b]\n');
        
        % This loop distribute the mapping value to each appropriate ROI
        for iBeta = 1:size(Betas,1)
            for iROI=1:size(ROIs,1)
                
                TEMP = AllMapping(ROIs{iROI,hs+2},:,iBeta);
                
                % This 2 following line are there to remove form further
                % analysis any vertex that would have at least one 0 or NaN
                % along its normal.
                TEMP(TEMP==0)=NaN;
                TEMP(any(isnan(TEMP),2),:)=nan(sum(any(isnan(TEMP),2)),NbLayers(iLayer));
                
                TEMP = TEMP';
                TEMP = TEMP(:);
                
                if all(isnan(TEMP))
                    TEMP = ones(size(TEMP))*Inf;
                end
                
                FeaturesAll{iROI,hs,iLayer}(strcmp(BetaList,Betas(iBeta)),:) = ...
                    repmat(TEMP', sum(strcmp(BetaList,Betas(iBeta))), 1);
                
                clear TEMP
                
            end
            
        end
        
        clear AllMapping
        
    end
    
end

if str2double(MatlabVer(1:4))>2013
if KillGcpOnExit
    delete(gcp)
end
else
end


fprintf('  Saving features matrices\n')

save(fullfile(AnalysisFolder, 'AllFeatures.mat'), 'FeaturesAll', '-v7.3');

for iLayer = 1:numel(NbLayers)
    
    for hs = 1:2
        
        if hs==1
            NameHS = '_lhs_';
        else
            NameHS = '_rhs_';
        end
        
        for iROI=1:size(ROIs,1)
            
            if ~all(FeaturesAll{iROI,hs,iLayer}(:)==Inf)
                
                if ~isempty(ROIs{iROI,hs+2}) && any(all(isnan(FeaturesAll{iROI,hs,iLayer}),2))
                    warning('Data for beta image %i is missing.\n', ...
                        str2double(unique(BetaList(all(isnan(FeaturesAll{iROI,hs,iLayer}),2)))))
                else
                    Features = FeaturesAll{iROI,hs,iLayer}; %#ok<NASGU>
                    save(fullfile(AnalysisFolder, ['Features_' ROIs{iROI,1} ...
                        NameHS num2str(NbLayers(iLayer)) '_Surf.mat']), 'Features', '-v7.3');
                end
                
            end
            
        end
    end
end

end