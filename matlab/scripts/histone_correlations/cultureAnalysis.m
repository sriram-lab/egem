function [solution] = cultureAnalysis(exp)
    %addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\metabolic_sensitivity')
    %addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\visualizations')
    
    addpath('/home/scampit/Desktop/egem/matlab/scripts/nutrient_sensitivity')
    addpath('/home/scampit/Desktop/egem/matlab/scripts/visualizations')
    
    load ./../../metabolic_models/eGEM_mm.mat
    load ./../../vars/ccle_geneExpression_vars.mat
    load ./../../vars/CCLE_Proteomics
    load allVars.mat

    reactions_of_interest = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
    model = eGEM_mm;
    BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
    model.c(BIOMASS_OBJ_POS) = 1;
    unique_cultures = unique(cultures);

    for cult = 1:length(unique_cultures)
        disp(unique_cultures(cult))
        culturePositions = find(ismember(string(cultures), string(unique_cultures(cult))));
        culturesCellNames = cell_names(culturePositions, :);

        cultureToGeneMatch = find(ismember(string(celllinenames_ccle1), ...
            string(culturesCellNames)));
        geneExpMatch = celllinenames_ccle1(cultureToGeneMatch);

        matched = intersect(string(culturesCellNames), string(geneExpMatch));
        matchedCulture = find(ismember(string(cell_names), string(matched)));
        matchedGenes = find(ismember(string(celllinenames_ccle1), string(matched)));
        culturesCellNames = cell_names(matchedCulture, :);
        geneExpMatch = celllinenames_ccle1(matchedGenes);

        cultureProteomicsValues = proteomics(matchedCulture, :);
        cultureMedium = medium(matchedCulture, :);

        [diffExp_genes] = find_diffexp_genes(model, geneExpMatch);

        for match = 1:length(culturesCellNames)

            ON_fieldname = string(strcat('ON_', geneExpMatch(match)));
            OFF_fieldname = string(strcat('OFF_', geneExpMatch(match)));
            str = strcat(string("ObjCoef = epsilon2_"), string(lower(cultureMedium(match))), string("(:, 1)"));
            eval(str);

            constrained_model = medium_LB_constraints(model, cultureMedium(match));

            reaction_name = char(reactions_of_interest(:, 1));
            reactions_to_optimize = find(ismember(constrained_model.rxns,...
                reaction_name));
            final_culture = string(unique_cultures(cult));
            
            switch exp
                case {'SRA', 'NoComp'}
                    kappa = 1;
                    rho = 1;
                    epsilon = 1E-3;
                    mode = 0; 
                    epsilon2 = 1E-3;
                    minfluxflag = true;

                    BIOMASS_OBJ_POS = find(ismember(constrained_model.rxns, 'biomass_objective')); 
                    constrained_model.c(BIOMASS_OBJ_POS) = 1;

                    for rxn = 1:length(reactions_to_optimize)
                        optimized_rxn = reactions_to_optimize(rxn);
                        constrained_model.c(optimized_rxn) = ObjCoef(rxn);

                        [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);
                        
                        fluxVarName = string(strcat(final_culture, '_flux_values(match, rxn) = solution.flux(optimized_rxn);'));
                        grateVarName = string(strcat(final_culture, '_grates(match, rxn) = solution.x(BIOMASS_OBJ_POS);'));
                        eval(fluxVarName);
                        eval(grateVarName);
                        constrained_model.c(reactions_to_optimize) = 0;
                        
                    end
                    
                case 'Comp'
                    kappa = 1;
                    rho = 1;
                    epsilon = 1E-3;
                    mode = 0; 
                    epsilon2 = 1E-3;
                    minfluxflag = true;

                    BIOMASS_OBJ_POS = find(ismember(constrained_model.rxns, 'biomass_objective')); 
                    constrained_model.c(BIOMASS_OBJ_POS) = 1;

                    reaction_positions = find(ismember(constrained_model.rxns, reactions_of_interest));
                    constrained_model.c(reaction_positions) = ObjCoef;

                    [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);

                    fluxVarName = string(strcat(final_culture, '_flux_values(match, :) = solution.flux(reaction_positions);'));
                    grateVarName = string(strcat(final_culture, '_grates(match, 1) = solution.x(BIOMASS_OBJ_POS);'));
                    eval(fluxVarName);
                    eval(grateVarName);
                    
                case 'FVA'
                    BIOMASS_OBJ_POS = find(ismember(constrained_model.rxns, 'biomass_objective')); 
                    constrained_model.c(BIOMASS_OBJ_POS) = 1;
                    reaction_positions = find(ismember(constrained_model.rxns, reactions_of_interest));
                    constrained_model.c(reaction_positions) = ObjCoef;
                    [~, maxFlux] = fluxVariability(constrained_model, 100, ...
                            'max', reactions_of_interest);
                    fluxVarName = string(strcat(final_culture, '_flux_values(match, :) = maxFlux;'));
                    eval(fluxVarName);                    
                    
            end
        end
        
        corrVarName = string(strcat(final_culture, '_corr = CultureCorr(', ...
            final_culture, '_flux_values, cultureProteomicsValues, final_culture, marks)'));
        eval(corrVarName)

        plotVarName = string(strcat('plotHistoneCorrelation(', final_culture, ...
            '_corr, "correlation", exp, "culture_corr")'));
        eval(plotVarName)

        histVarName = strcat("makeHist(", final_culture, ...
            "_flux_values, string(unique_cultures(cult)), 'culture_hist', exp);");
        eval(histVarName)
    end
end

 