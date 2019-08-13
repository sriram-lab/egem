% testDMLeaks finds reactions that contain metabolite leaks for demand
% reactions.
% @author: Scott Campit
function [DMChecks] = testDMLeaks(model)
    
    [selectedExchangeReactions] = getAllExchangeRxns(model);
    DemandReactions = strmatch('DM_', selectedExchangeReactions);

    [LeakRxnsDM, ~, LeakRxnsFluxVectorDM] = ...
        fastLeakTest(model, model.rxns(DemandReactions), 'true');
    
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