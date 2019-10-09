% testEXLeaks finds reactions that contain metabolite leaks for exchanged
% reactions.
% @author: Scott Campit
function [EXChecks] = testEXLeaks(model)
    [selectedExchangeReactions] = getAllExchangeRxns(model);
    ExchangeReactions = strmatch('EX_', selectedExchangeReactions);

    [LeakEXRxns, ~, LeakEXRxnsFluxVector] = ...
        fastLeakTest(model, model.rxns(ExchangeReactions),'false');
    
    disp(LeakEXRxns)
    
    if ~isempty(LeakEXRxns)
        warning('Model leaks metabolites!')
    else
        disp('Model has no metabolite leaks!')
    end
    
    fields = { ...
        'Leak_EX_Rxns', 'LeakEXRxnsFluxVector' ...
    };
    values = { ...
        LeakEXRxns, LeakEXRxnsFluxVector ...
    };
    
    for i = 1:length(fields)
        EXChecks.(fields{i}) = values{i};
    end
    
end