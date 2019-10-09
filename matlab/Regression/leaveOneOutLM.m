function leaveOneOutLM(X, Y)
    for feat = 1:size(X, 2)
        tmpArr = X;
        tmpArr(:, feat) = [];
        tmpModel = fitlm(tmpArr, Y);
        tmpModel
        anova(tmpModel)
        
%         figure(1);
%         plotDiagnostics(tmpModel, 'cookd')
%         
%         figure(2);
%         plotResiduals(tmpModel, 'fitted')
%         
%         figure(3);
%         plotResiduals(tmpModel, 'probability')
        
    end
end