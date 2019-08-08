% `calc_metabolic_metrics`: Compute the solution from flux balance analysis
% and flux variability analysis. 
% @author: Scott Campit

function [solution] = calc_metabolic_metrics(model, reaction_positions, ...
    metabolites_of_interest, grate, sense, fva_rxns, objCoef, exp)

switch exp
    case {'sra', 'no_competition', 'competition'}
        model.c(reaction_positions) = objCoef;
        gurobi_soln = optimizeCbModel(model);
        solution.reactionFlux = gurobi_soln.v(reaction_positions);
        solution.reducedCosts = gurobi_soln.w(reaction_positions);
        solution.shadowPrice = gurobi_soln.y(find(ismember(model.mets, ...
            metabolites_of_interest)));
        solution.maxFlux = [];
        solution.minFlux = [];
        
        if reaction_positions == find(ismember(model.rxns, ...
                'biomass_objective'))
            model.c(reaction_positions) = 1;
        else
            model.c(reaction_positions) = 0;
        end

        
        
    case 'fva'
        [minFlux, maxFlux] = fluxVariability(model, grate, sense, fva_rxns);
        solution.name = 'FluxVariabilityAnalysis';
        solution.minflux = minFlux;
        solution.maxflux = maxFlux;
        solution.reactionFlux = [];
        solution.shadowPrice = [];
        solution.reducedCosts = [];
    
end

end