function v = testOptimization(model, listOfReactions)
    for i = 1:length(listOfReactions)
        rxnOfInterest = find(ismember(model.rxns, listOfReactions(i)));
        model.c(rxnOfInterest) = 1E-3;
        soln = optimizeCbModel(model);
        v(i) = soln.x(rxnOfInterest);
        model.c(rxnOfInterest) = 0;
    end
end