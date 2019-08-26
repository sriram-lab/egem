function [deadEndMets, deadEndRxns, RxnFormula] = findDeadEnds(model)
    deadEndPos = detectDeadEnds(model);
    deadEndMets = model.mets(deadEndPos);
    
    [deadEndRxns, RxnFormula] = findRxnsFromMets(model, deadEndMets);
    
end