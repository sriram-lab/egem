function [flux_values, grates] = singleReactionAnalysis(model, objCoef)
    
    kappa = 1;
    rho = 1;
    epsilon = 1E-3;
    mode = 0; 
    epsilon2 = 1E-3;
    minfluxflag = true;
    
    BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
    model.c(BIOMASS_OBJ_POS) = 1;
    reactions_to_optimize = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
    
    for rxn = 1:length(reactions_to_optimize)
        optimized_rxn = reactions_to_optimize(rxn);
        constrained_model.c(reactions_to_optimize) = objCoef(rxn);

        [~, solution] =  constrain_flux_regulation...
           (constrained_model, diffExp_genes.(ON_fieldname), ...
            diffExp_genes.(OFF_fieldname), ...
            kappa, rho, epsilon, mode, [], minfluxflag);

        flux_values(match, rxn) = solution.flux(optimized_rxn);
        grates(match, rxn) = solution.flux(BIOMASS_OBJ_POS);
        constrained_model.c(reactions_to_optimize) = 0;
    end
    
end