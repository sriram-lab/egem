function linear_optimization_experiments(reactions_of_interest)

    function single_reaction_analysis(reactions_of_interest)
        for reaction = 1:length(reactions_of_interest)
            [solution] = calc_metabolic_metrics(excess_model,...
                            objectivePosition, metabolites_of_interest,...
                            [], [], [], 1, exp);
            
            
            