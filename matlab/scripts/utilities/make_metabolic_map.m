function [map] = make_metabolic_map(metabolicmodel)
    rxnID = metabolicmodel.rxns;
    rxnName = metabolicmodel.rxnNames;
    rxnLB = num2cell(metabolicmodel.lb);
    rxnUB = num2cell(metabolicmodel.ub);
    map = [rxnID, rxnName, rxnLB, rxnUB];
    map(1,:) = {'BiGG ID', 'Reaction Name', 'LB', 'UB'};
end