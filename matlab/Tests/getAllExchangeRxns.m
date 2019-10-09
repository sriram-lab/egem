% @author: Scott Campit
function [UniqueExchangeReactions] = getAllExchangeRxns(model)
    
    testModel = model;
    
    ExchangeReactions = strmatch('EX_', testModel.rxns);
    DemandReactions = strmatch('DM_', testModel.rxns);
    
    AllExchangeRxns = (find(full((sum(abs(testModel.S) == 1, 1) == 1) & ...
        (sum(testModel.S ~= 0) == 1))))';
    
    UniqueExchangeReactions = unique([ExchangeReactions; ...
        DemandReactions; AllExchangeRxns]);
end