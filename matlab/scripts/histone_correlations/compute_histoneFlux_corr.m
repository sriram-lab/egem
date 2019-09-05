% `compute_histoneFlux_corr` computes the correlation between histone
% proteomics data and metabolic fluxes computed by the eGEM
% @author: Scott Campit
function [solution] = compute_histoneFlux_corr(eGEM, reactions_of_interest, ...
    eps_struct, mode, epsilon, rho, kappa, minfluxflag, exp, ...
    proteomics_dataset, fva_grate)

    [cell_names, medium, marks, proteomics] = load_proteomics('CCLE');
        
    load ./../vars/supplementary_software_code...
        celllinenames_ccle1... % CCLE cellline names
        ccleids_met... % Gene symbols
        ccle_expression_metz % Z-transformed gene expression

    proteomics_CL_match = find(ismember(cell_names, celllinenames_ccle1));
    BIOMASS_OBJ_POS = find(ismember(eGEM.rxns, 'biomass_objective')); 

    proteomics = proteomics(proteomics_CL_match, :);
    cell_names = cell_names(proteomics_CL_match,1);
    geneExp_CL_match = find(ismember(celllinenames_ccle1, cell_names));

    [diffExp_genes] = find_diffexp_genes(eGEM, geneExp_CL_match);

    number_of_matched_proteomics = length(proteomics_CL_match);
    for match = 1:number_of_matched_proteomics

        constrained_model = eGEM;
        constrained_model = media(constrained_model, medium(match));
        obj_coef = 1E-3; % This is temporary. I will have to recalculate these values based on the new formulation.
        constrained_model.c(BIOMASS_OBJ_POS) = 1;

        reaction_name = char(reactions_of_interest(:, 1));
        switch exp
            case 'non-competitive_cfr'
                for rxn = 1:length(reactions_of_interest(:,1))
                    model_to_test = constrained_model;
                    rxnpos = [find(ismember(model_to_test.rxns, reactions_of_interest(rxn)))];
                    model_to_test.c(rxnpos) = obj_coef(rxn, 1);
                    [flux, ~, ~] = constrain_flux_regulation(model_to_test,  ...
                        onreactions, offreactions, kappa, rho, epsilon, mode, [], ...
                        minfluxflag);
                    all_flux_values(match, rxn) = flux(rxnpos);
                    model_to_test.c(rxnpos) = 0;
                end
            case 'competitive_cfr'
                model_to_test = constrained_model;
                rxnpos = [find(ismember(model_to_test.rxns, reaction_name))];
                model_to_test.c(rxnpos) = obj_coef(:, 1);
                [flux, ~, ~] =  constrain_flux_regulation(model_to_test,...
                    onreactions, offreactions, kappa, rho, epsilon, mode , [], ...
                    minfluxflag);
                all_flux_values(match,:) = flux(rxnpos);
            case 'fva'
                model_to_test = constrained_model;
                rxnpos = [find(ismember(model_to_test.rxns, reaction_name))];
                model_to_test.c(rxnpos) = obj_coef(:, 1);
                [~, ~, ~, ~, flux, ~] =...
                    calc_metabolic_metrics(model_to_test, [], [], fva_grate,...
                    'max', reactions_of_interest, [], 'fva');
                all_flux_values(match,:) = flux;
            case 'tissueAnalysis'
                unique_tissues = unique(tissues);
                for tiss = 1:length(unique_tissues)
                    tissue = unique_tissues(tiss);
                    tissue_positions = find(ismember(tissues, tissue));
                    tissue_proteomics = proteomics(tissue_psoitions, :);

            case 'mediumAnalysis'
        end
    end

    [rho, pval] = corr(all_flux_values, proteomics);
    rxns = reactions_of_interest(:, 3);

    %% Save data in struct
    solution = struct('Name', proteomics_dataset);

    fields = {...
            'HistoneMark'; 'Reaction'; ...
            'R'; 'Pvalue'; 'Flux'; 'Proteomics'
        };

    values = {...
        marks; rxns; ...
        rho; pval; all_flux_values; proteomics
        };

    for match=1:length(fields)
        solution.(fields{match}) = values{match};
    end

end