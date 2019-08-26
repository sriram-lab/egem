% testDMLeaks finds reactions that contain metabolite leaks for demand
% reactions.
% @author: Scott Campit
function [DMChecks] = testDMLeaks(model)
    
    [selectedExchangeReactions] = getAllExchangeRxns(model);
    allExchangeReactions = model.rxns(selectedExchangeReactions);
    DemandReactions = allExchangeReactions(strmatch('DM_', allExchangeReactions));
    DM_pos = find(ismember(model.rxns, DemandReactions));

    [LeakRxnsDM, ~, LeakRxnsFluxVectorDM] = ...
        fastLeakTest(model, model.rxns(DM_pos), 'false');
    
    disp(LeakRxnsDM)
    
    if ~isempty(LeakRxnsDM)
        disp('SUCCESS: Leak free when demand reactions are added!')
    else
         warning('WARNING: Model contains leaky Demand Reactions')
    end
    
    fields = { ...
        'Leak_DM_Rxns', 'LeakDMRxnsFluxVector' ...
    };
    values = { ...
        LeakRxnsDM, LeakRxnsFluxVectorDM ...
    };
    
    for i = 1:length(fields)
        DMChecks.(fields{i}) = values{i};
    end
end