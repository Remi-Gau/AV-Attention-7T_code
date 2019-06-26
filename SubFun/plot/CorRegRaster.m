function [rho,slope]=CorRegRaster(Profiles, DesMat, iToPlot, X_sort)
% computes vertex wise correlation / regression between predictor and predicted rasters
% plots


for iPerc = 1:size(Profiles,1)
    Y = squeeze(Profiles(iPerc,:,:));
    
    X = [];
    for iSubj=1:size(Y,2)
        X((1:6)+6*(iSubj-1),(1:size(DesMat,2))+size(DesMat,2)*(iSubj-1)) = DesMat;
    end
    
    Y = Y(:);
    B = pinv(X)*Y;
    
    Cst_tmp = B(1:size(DesMat,2):size(X,2),:);
    Lin_tmp = B(2:size(DesMat,2):size(X,2),:);
    
    if iToPlot==1
        Y_sort(:,iPerc,:)=Cst_tmp;
    else
        Y_sort(:,iPerc,:)=Lin_tmp;
    end
    
end

for iSubj=1:size(Y_sort,1)
    R=corrcoef(X_sort(iSubj,:),Y_sort(iSubj,:));
    rho(iSubj) = R(1,2);
    beta = glmfit(X_sort(iSubj,:), Y_sort(iSubj,:), 'normal');
    slope(iSubj) = beta(2);
end
end