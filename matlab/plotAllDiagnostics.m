function plotAllDiagnostics(model)
    disp("ANOVA Results")
    anova(model)
    
    plotDiagnostics(model, 'cookd')
    
    plotResiduals(model)
    
    plotResiduals(model, 'fitted')
end
