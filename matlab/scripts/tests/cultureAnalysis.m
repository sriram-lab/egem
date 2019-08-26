%initCobraToolbox
%changeCobraSolver('gurobi');
function cultureAnalysis()
    %addpath('/home/scampit/Desktop/eGEM/matlab/scripts/metabolic_sensitivity/')
    addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\metabolic_sensitivity')

    %addpath('/home/scampit/Desktop/eGEM/matlab/scripts/visualizations')
    addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\visualizations')
    addpath('C:\Users\scampit\Desktop\eGEM\matlab\scripts\metabolic_sensitivity')
    load ./../../metabolic_models/eGEM_mm.mat
    load ./../../vars/ccle_geneExpression_vars.mat
    load ./../../vars/CCLE_Proteomics

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
            obj_coef = [1E-3, 1E-3, 1E-3, 1E-3];

            ON_fieldname = string(strcat('ON_', geneExpMatch(match)));
            OFF_fieldname = string(strcat('OFF_', geneExpMatch(match)));

            constrained_model = medium_LB_constraints(model, cultureMedium(match));

            reaction_name = char(reactions_of_interest(:, 1));
            reactions_to_optimize = [find(ismember(constrained_model.rxns,...
                reaction_name))];

            for rxn = 1:length(reactions_to_optimize)
                optimized_rxn = reactions_to_optimize(rxn);
                constrained_model.c(reactions_to_optimize) = obj_coef(rxn);

                kappa = 1;
                rho = 1;
                epsilon = 1E-3;
                mode = 0; 
                epsilon2 = 1E-3;
                minfluxflag = true;

                [cellLine_model, solution] =  constrain_flux_regulation...
                   (constrained_model, diffExp_genes.(ON_fieldname), ...
                    diffExp_genes.(OFF_fieldname), ...
                    kappa, rho, epsilon, mode, [], minfluxflag);

                medium_flux_values(match, rxn) = solution.flux(optimized_rxn);
                medium_grates(match, rxn) = solution.flux(BIOMASS_OBJ_POS);
                constrained_model.c(reactions_to_optimize) = 0;
                final_culture = string(unique_cultures(cult));
                fluxVarName = string(strcat(final_culture, '_flux_values(match, rxn) = solution.flux(optimized_rxn);'));
                grateVarName = string(strcat(final_culture, '_grates(match, rxn) = solution.x(BIOMASS_OBJ_POS);'));
                eval(fluxVarName);
                eval(grateVarName);
            end
        end

        corrVarName = string(strcat(final_culture, '_corr = ', 'CultureCorr(', ...
            final_culture, '_flux_values, mediaProteomicsValues, mediaList(med), marks)'));
        eval(corrVarName)

        plotVarName = string(strcat('plot_heatmap(', final_culture, ...
            '_corr', ', [], "correlation", [], [], [])'));
        eval(plotVarName)
end

 