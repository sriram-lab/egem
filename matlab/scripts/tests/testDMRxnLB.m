function [wrongLB] = testDMRxnLB(model)
    DMlb = find(model.lb(strmatch('DM_', model.rxns)) < 0);
    
    if isempty(DMlb)
        disp('SUCCESS: No demand reaction can have flux in backward direction.')
    else
        warning('Demand reaction can have flux in backward direction.')
    end
    
    wrong_LB = model.rxns(strmatch('DM_', model.rxns))
    
    fields = {'wrongLB'};
    values = {wrong_LB};
    
    for i = 1:length(fields)
        wrongLB.(fields{i}) = values{i};
    end
end
