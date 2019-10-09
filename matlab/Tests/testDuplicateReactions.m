% @author: Scott Campit
function [DuplicateChecks] = testDuplicateReactions(model)
    method = 'S';
    removeFlag = 0;
    [~, removedRxnInd, ~] = checkDuplicateRxn(model, ...
        method, removeFlag, 0);
    
    duplicateRxn = model.rxns(removedRxnInd);
    actual_duplicateRxn = ~strmatch('copy', duplicateRxn);
    
    if isempty(actual_duplicateRxn)
        disp('SUCCESS: No duplicated reactions in model.')
    else
        warning('WARNING: Duplicated reactions in model.')
    end
    
    fields = {'DuplicatedReactions'};
    values = {actual_duplicateRxn};
    
    for i = 1:length(fields)
        DuplicateChecks.(fields{i}) = values{i};
    end
    
end