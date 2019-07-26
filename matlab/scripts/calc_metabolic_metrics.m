%% Get metabolic model outputs
function [soln, flux, reducedCosts, shadowPrice, maxFlux, minFlux] = ...
    calc_metabolic_metrics(input_model, rxnPos, metabolites_of_interest, ...
    grate, sense, fva_rxns, objCoef, exp)

switch exp
    case {'sra', 'no_competition', 'competition'}
        input_model.c(rxnPos) = objCoef;
        soln = optimizeCbModel(input_model);
        flux = soln.v(rxnPos);
        reducedCosts = soln.w(rxnPos);
        shadowPrice = soln.y(find(ismember(input_model.mets, metabolites_of_interest)));
        
        if rxnPos == find(ismember(input_model.rxns, 'biomass_objective'))
            input_model.c(rxnPos) = 1;
        else
            input_model.c(rxnPos) = 0;
        end
        
        maxFlux = 0;
        minFlux = 0;

    case 'fva'
        [minFlux, maxFlux] = fluxVariability(input_model, grate, sense, fva_rxns);
        soln = 0;
        flux = 0;
        shadowPrice = 0;
        reducedCosts = 0;
end

end