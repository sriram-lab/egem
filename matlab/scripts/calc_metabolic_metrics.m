%% Get metabolic model outputs
function [soln, v, rc, sp, maxFlux, minFlux] = calc_metabolic_metrics(model,...
        reactions_of_interest, metabolites_of_interest, grate, sense,...
        fva_rxns, objCoef, type)
switch type
    case {'sra', 'no_competition'}
        model.c(reactions_of_interest) = objCoef;
        soln = optimizeCbModel(model);
        v = soln.v(reactions_of_interest);
        rc = soln.w(reactions_of_interest);
        sp = soln.y(find(ismember(model.mets, metabolites_of_interest)));
        if reactions_of_interest == find(ismember(model.rxns, 'biomass_objective'))
          model.c(reactions_of_interest) = 1;
        else
          model.c(reactions_of_interest) = 0;
        end
        maxFlux=0;
        minFlux=0;

    case 'competition'
        model.c(reactions_of_interest) = objCoef;
        soln = optimizeCbModel(model);
        v = soln.v(reactions_of_interest);
        rc = soln.w(reactions_of_interest);
        sp = soln.y(find(ismember(model.mets, metabolites_of_interest)));
        maxFlux=0;
        minFlux=0;

    case 'fva'
        model.c(reactions_of_interest) = objCoef;
        [minFlux, maxFlux] = fluxVariability(model, grate, sense,...
            fva_rxns);
        v=0;
        sp=0;
        rc=0;
        soln=0;

end
end
