function Beta = PlotScatterDensity2(X,Y,XRange,YRange,Bins)
    
%     sum(isnan(X))
%     sum(isnan(Y))

    Xi = linspace(XRange(1),XRange(2),Bins);
    Yi = linspace(YRange(1),YRange(2),Bins);
        
    Xr = interp1(Xi, 1:numel(Xi), X, 'nearest'); clear Xi
    Yr = interp1(Yi, 1:numel(Yi), Y, 'nearest'); clear Yi
    
    sum(isnan(Xr));
    sum(isnan(Yr));

    Yr(isnan(Xr))=[];
    Xr(isnan(Xr))=[];
    
    Xr(isnan(Yr))=[];
    Yr(isnan(Yr))=[];
    
    PDF = accumarray([round(Yr) ; round(Xr)]', 1, [Bins Bins]);
    % PDF=PDF/max(PDF(:));
    
    B = glmfit(Xr, Yr, 'normal');
    
    A=colormap;
    A(1,:) = [1 1 1];
    colormap(A);
    imagesc(PDF);
    axis square

    plot([0 Bins], ([0 Bins])*B(2)+B(1), 'k', 'LineWidth', 2);
    plot([Bins/2 Bins/2], [0 Bins], 'k')
    plot([0 Bins], [Bins/2 Bins/2], 'k')
    
    Y(isnan(X))=[]; X(isnan(X))=[];
    X(isnan(Y))=[]; Y(isnan(Y))=[];
    
    Beta = glmfit(X, Y, 'normal');
    
    R=corrcoef(X,Y);
    Beta(4) = R(1,2);
    
    X = (X-mean(X))/std(X);
    Y = (Y-mean(Y))/std(Y);

    R=corrcoef(X,Y);
    Beta(3) = R(1,2);
    
    if abs(Beta(4)-Beta(3))>0.01
        fprintf('rho=%3.3f ; rho_Z=%3.3f', Beta(4), Beta(3))
    end
    
    t = text(10, 20, sprintf('B_1=%3.3f; B_0=%3.3f; rho=%3.3f', Beta(2), Beta(1), Beta(3)));
    set(t,'fontsize',14)
    
    set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'fontsize', 6, ...
        'xtick', linspace(1, Bins, 6) ,'xticklabel', linspace(XRange(1),XRange(2), 6),...
        'ytick', linspace(1, Bins, 6) ,'yticklabel', linspace(YRange(1),YRange(2), 6))
    
    axis([0 Bins 0 Bins])
    
    box off

end

