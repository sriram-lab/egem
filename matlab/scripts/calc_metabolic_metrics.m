%% Get metabolic model outputs 
function [v, rc, sp, maxFlux, minFlux] = calc_metabolic_metrics(model,...
        reactions_of_interest, metabolites_of_interest, grate, sense,...
        fva_rxns, objcoef, type)
switch type
    case 'sra'
        model.c(reactions_of_interest) = objcoef;
        soln = optimizeCbModel(model);
        v = soln.v(reactions_of_interest);
        rc = soln.w(reactions_of_interest);
        sp = soln.y(find(ismember(model.mets, metabolites_of_interest)));
        model.c(reactions_of_interest) = 0;
        maxFlux=0;
        minFlux=0;
        
    case 'fba'
        model.c(reactions_of_interest) = objcoef;
        soln = optimizeCbModel(model);
        v = soln.v(reactions_of_interest);
        rc = soln.w(reactions_of_interest);
        sp = soln.y(find(ismember(model.mets, metabolites_of_interest)));
        %model.c(reactions_of_interest) = 0;
        maxFlux=0;
        minFlux=0;
        
    case 'fva'
        model.c(reactions_of_interest) = objcoef;
        [minFlux, maxFlux] = fluxVariability(model, grate, sense,...
            fva_rxns);
        %model.c(reactions_of_interest) = 0;
        v=0;
        sp=0;
        rc=0;
        
end
end