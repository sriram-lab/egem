function [grate, flux, solution] = competitiveFBA(model, objCoef)
    
    BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
    model.c(BIOMASS_OBJ_POS) = 1;
    
    reactions_to_optimize = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
    reaction_positions = find(ismember(model.rxns, reactions_to_optimize));
    model.c(reaction_positions) = objCoef;
    
    [~, solution] =  constrain_flux_regulation...
           (constrained_model, diffExp_genes.(ON_fieldname), ...
            diffExp_genes.(OFF_fieldname), ...
            kappa, rho, epsilon, mode, [], minfluxflag);
    
    grate = FBAsoln.x(BIOMASS_OBJ_POS);
    flux = FBAsoln.x(reaction_positions);
    
end