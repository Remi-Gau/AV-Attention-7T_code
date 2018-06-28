
function [Betas, FilesList] = CorresBetaVTK(Folder, CopyFiles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
%     Folder = '/data/AV_Integration_2/Subjects_Data/Subject_06/Transfer/BetaMapping/exp-0000';
    SubjectList = ['07'];
    for i=1:1;
        CorresBetaVTK(fullfile(pwd, 'Subjects_Data', ['Subject_' SubjectList(i,:)],...
            'Transfer', 'BetaMapping','exp-0000'), 1)
    end
end

if nargin<2
    CopyFiles = 1;
end

StartFolder = pwd;

cd(Folder)

FolderList = dir;

IsDir = find([FolderList.isdir]);


% Betas = cell(length(IsDir)-2, 1);
% FilesList = cell(length(IsDir)-2, 2, 1);
Betas = {};
FilesList = {};


irow = 1;

for iFolder = 3:length(IsDir)
    
    cd(fullfile(Folder, FolderList(IsDir(iFolder)).name))
    
    temp=dir;
    
    if strcmp(temp(find([temp.isdir])>2).name, 'SurfaceMeshMapping')
        
        
        FileList = dir('*.input');
        
        Side = NaN;
        Surf = NaN;       
        
        FileContent = fileread(FileList.name);
        
        
        if isempty(strfind(FileContent, 'rbeta_'))
            error('Mapping was not done on a beta image!!')
        end
        
        if ~isempty(strfind(FileContent, 'lcr_gm_avg'))
            Side=1; Surf=1; SideInsert='_lcr';
        elseif ~isempty(strfind(FileContent, 'rcr_gm_avg'))
            Side=2; Surf=1; SideInsert='_rcr';
        end
        
        if ~isnan(Side) && Surf==1
            
            B = strfind(FileContent, 'rbeta_');
            
            B = FileContent(B(1)+6:B(1)+9);
            
            File = dir(fullfile(Folder, FolderList(IsDir(iFolder)).name, ...
               'SurfaceMeshMapping', '*.vtk'));
            
            if CopyFiles
                copyfile(fullfile(Folder, FolderList(IsDir(iFolder)).name, ...
                'SurfaceMeshMapping', File.name), ...
                fullfile(Folder, ['Beta_' B SideInsert '.vtk']));
            end

        end

        irow = irow + 1;
        
        cd ..
        
    end
    
end

cd(StartFolder)

end