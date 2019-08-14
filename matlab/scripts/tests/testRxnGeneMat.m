function [emptyRxnColumnPosition] = testRxnGeneMat(model)
    E = find(all(model.rxnGeneMat == 0));
    if isempty(E)
        disp('SUCCESS: No empty columns in model.rxnGeneMat')
    else
        warning('WARNING: Empty columns in model.rxnGeneMat')
    end
    
    emptyRxns = model.rxns(E)
    
    fields = {'EmptyRxns'};
    values = {emptyRxns};
    
    for i = 1:length(fields)
        emptyRxnColumnPosition.(fields{i}) = values{i};
    end
end