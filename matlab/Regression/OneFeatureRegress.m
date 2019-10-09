function OneFeatureRegress(X, Y)
    i = 1;
    for col = 1:size(X, 2)
        tmp = X(:, col);
        tmpModel = fitlm(tmp, Y);
        tmpModel
        anova(tmpModel)
        
        %figure(i);
        %plotDiagnostics(tmpModel, 'cookd')
        %i = i+1;
        
        %figure(i);
        %plotResiduals(tmpModel, 'fitted')
        %i = i+1;
        
        %figure(i);
        %plotResiduals(tmpModel, 'probability')
        %i = i+1;

    end
end