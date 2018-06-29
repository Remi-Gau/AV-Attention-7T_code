clc; clear; close all;

SubjectsList = [...
    '02'
    '03'
    '04'
    '06'
    '07'
    '08'
    '09'
    '11'
    '12'
    '13'
    '14'
    '15'
    '16'];
NbSubj = size(SubjectsList,1);

StartFolder = pwd;

NbBins = 128;

HistT1 = zeros(NbBins,NbSubj);
HistUNI = zeros(NbBins,NbSubj);
HistINV2 = zeros(NbBins,NbSubj);

RangeT1 = zeros(NbSubj,2);
RangeUNI = zeros(NbSubj,2);
RangeINV2 = zeros(NbSubj,2);

for iSubj=1:NbSubj
    
    fprintf('\n\nSubject %s \n', SubjectsList(iSubj,:))
    
    SubjFolder = fullfile(StartFolder, 'Subjects_Data', ['Subject_' SubjectsList(iSubj,:)]);
    
    cd(SubjFolder)
    
    T1 = spm_read_vols(spm_vol('T1.nii'));
    HistT1(:,iSubj) = histc(T1(:),linspace(0,4000,NbBins));
    RangeT1(iSubj,:) = [min(T1(:)) max(T1(:))];
    
    UNI = spm_read_vols(spm_vol('UNI.nii'));
    HistUNI(:,iSubj) = histc(UNI(:),linspace(0,4100,NbBins));
    RangeUNI(iSubj,:) = [min(UNI(:)) max(UNI(:))];
    
    INV2 = spm_read_vols(spm_vol('INV2.nii'));
    HistINV2(:,iSubj) = histc(INV2(:),linspace(0,3500,NbBins));
    RangeINV2(iSubj,:) = [min(INV2(:)) max(INV2(:))];
    
end

RangeT1
RangeUNI
RangeINV2

cd(StartFolder)

save('HistData.mat', 'SubjectsList', 'HistT1', 'HistUNI', 'HistINV2', ...
    'RangeT1', 'RangeUNI', 'RangeINV2')


%%
SubjFolder = '/Users/rxg243/Dropbox/PhD/Experiments/Magdeburg';

cd(SubjFolder)

NbBins = 128;

% T1 = spm_read_vols(spm_vol('T1.nii'));
T1 = spm_read_vols(spm_vol('T1_histMatch.nii'));
MagdeHistT1 = histc(T1(:),linspace(0,4000,NbBins));
MagdeRangeT1(1,:) = [min(T1(:)) max(T1(:))];

% UNI = spm_read_vols(spm_vol('UNI.nii'));
UNI = spm_read_vols(spm_vol('UNI_histMatch.nii'));
MagdeHistUNI = histc(UNI(:),linspace(0,4100,NbBins));
MagdeRangeUNI(1,:) = [min(UNI(:)) max(UNI(:))];

% INV2 = spm_read_vols(spm_vol('INV2.nii'));
INV2 = spm_read_vols(spm_vol('INV2_histMatch.nii'));
MagdeHistINV2 = histc(INV2(:),linspace(0,3500,NbBins));
MagdeRangeINV2(1,:) = [min(INV2(:)) max(INV2(:))];

MagdeRangeT1
MagdeRangeUNI
MagdeRangeINV2

% Color for Subjects
COLOR_Subject= [
    0,0,0;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    0,0,130;
    177,89,40;
    125,125,125];
COLOR_Subject=COLOR_Subject/255;

StartFolder = '/Users/rxg243/Dropbox/PhD/Experiments/7T_Project/Final_2';
cd(StartFolder)

load('HistData.mat', 'SubjectsList', 'HistT1', 'HistUNI', 'HistINV2')

NbSubj = size(SubjectsList,1);



figure(1)
subplot(231)
hold on
% for iSubj=1:NbSubj
%     plot(HistT1(:,iSubj), 'color', COLOR_Subject(iSubj,:))
% end
shadedErrorBar(1:NbBins,mean(HistT1,2),std(HistT1,[],2),'-b',1)
plot(MagdeHistT1, 'r')
set(gca,'YScale','log', 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,4000,5));
t = title('T1 map');
set(t, 'FontSize', 16);

subplot(232)
hold on
% for iSubj=1:NbSubj
%     plot(HistUNI(:,iSubj), 'color', COLOR_Subject(iSubj,:))
% end
shadedErrorBar(1:NbBins,mean(HistUNI,2),std(HistUNI,[],2),'-b',1)
plot(MagdeHistUNI, 'r')
set(gca,'YScale','log', 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,4100,5));
t = title('Weighted T1');
set(t, 'FontSize', 16);

subplot(233)
hold on
for iSubj=1:NbSubj
%     plot(HistINV2(:,iSubj), 'color', COLOR_Subject(iSubj,:))
    plot(HistINV2(:,iSubj), 'b')
end
% shadedErrorBar(1:NbBins,mean(HistINV2,2),std(HistINV2,[],2),'-b',0)
plot(MagdeHistINV2, 'r', 'linewidth', 2)
set(gca,'YScale','log', 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,3500,5));
t = title('INV2');
set(t, 'FontSize', 16);


subplot(234)
hold on
% for iSubj=1:NbSubj
%     plot(HistT1(:,iSubj), 'color', COLOR_Subject(iSubj,:))
% end
shadedErrorBar(1:NbBins,mean(HistT1,2),std(HistT1,[],2),'-b',0)
plot(MagdeHistT1, 'r')
set(gca, 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,4000,5));

subplot(235)
hold on
% for iSubj=1:NbSubj
%     plot(HistUNI(:,iSubj), 'color', COLOR_Subject(iSubj,:))
% end
shadedErrorBar(1:NbBins,mean(HistUNI,2),std(HistUNI,[],2),'-b',0)
plot(MagdeHistUNI, 'r')
set(gca, 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,4100,5));

subplot(236)
hold on
for iSubj=1:NbSubj
%     plot(HistINV2(:,iSubj), 'color', COLOR_Subject(iSubj,:))
    plot(HistINV2(:,iSubj), 'b')
end
% shadedErrorBar(1:NbBins,mean(HistINV2,2),std(HistINV2,[],2),'-b',0)
plot(MagdeHistINV2, 'r', 'linewidth', 2)
set(gca, 'XTick',linspace(0,NbBins,5),'XTickLabel',linspace(0,3500,5));

return

%%
MeanT1=mean(HistT1,2);
MeanUNI=mean(HistUNI,2);
MeanINV2=mean(HistINV2,2);


    
for iSubj=1:size(SubjectsList,1)
    
    ErrorT1(iSubj) = sum((HistT1(:,iSubj)-MeanT1).^2);
    ErrorUNI(iSubj) = sum((HistUNI(:,iSubj)-MeanUNI).^2);
    ErrorINV2(iSubj) = sum((HistINV2(:,iSubj)-MeanINV2).^2);
    
end


find(ErrorT1==min(ErrorT1))
find(ErrorUNI==min(ErrorUNI))
find(ErrorINV2==min(ErrorINV2))