initCobraToolbox
changeCobraSolver('gurobi');

load ./../models/eGEM.mat 
model = eGEM;

path = './../data/';
proteomics_array = readcell(strcat(path, 'GCP_proteomics_remapped.csv'));

cell_names = proteomics_array(2:end, 1);
tissues = proteomics_array(2:end, 2);
marks = proteomics_array(1, 3:44);
medium = proteomics_array(2:end, 45);

proteomics = proteomics_array(2:end, 3:44);
missing_values = cellfun(@ismissing, proteomics);
proteomics(missing_values) = {NaN};
proteomics = cell2mat(proteomics);
proteomics = knnimpute(proteomics);

load ./../vars/supplementary_software_code...
    celllinenames_ccle1... % CCLE cellline names
    ccleids_met... % Gene symbols
    ccle_expression_metz % Z-transformed gene expression

switch in_depth_analysis
    case 'tissueAnalysis'
        unique_tissues = unique(tissues);
        
        for tiss = 1:length(unique_tissues)
            BIOMASS_OBJ_POS = find(ismember(eGEM.rxns, 'biomass_objective')); 
            
            % match against each tissue
            tissue = unique_tissues(tiss);
            disp(tissue)
            tissue_positions = find(ismember(string(tissues), string(tissue)));
            tissue_proteomics = proteomics(tissue_positions, :);
            tissue_cellNames = cell_names(tissue_positions, :);
            tissue_medium = medium(tissue_positions, :);
            
            % match against the CCLE gene expression dataset
            proteomics_CL_match_positions = find(ismember(string(tissue_cellNames), ...
                string(celllinenames_ccle1)));

            tissue_matched_proteomics = tissue_proteomics(proteomics_CL_match_positions, :);
            tissue_matched_cellNames = tissue_cellNames(proteomics_CL_match_positions,1);
            geneExp_CL_match_positions = find(ismember(string(celllinenames_ccle1),...
                string(tissue_matched_cellNames)));
            geneExp_CL_match_cellNames = celllinenames_ccle1(geneExp_CL_match_positions);
            [diffExp_genes] = find_diffexp_genes(eGEM, geneExp_CL_match_cellNames);
        end
            for match = 1:length(proteomics_CL_match_positions)
                obj_coef = eps_struct.RPMI;
                constrained_model = medium_LB_constraints(eGEM, tissue_medium(match));
                constrained_model.c(BIOMASS_OBJ_POS) = 1;
                reaction_name = char(reactions_of_interest(:, 1));
                %constrained_model.lb(constrained_model.lb < -1000) = -1000;
                reactions_to_optimize = [find(ismember(constrained_model.rxns,...
                    reaction_name))];
                constrained_model.c(reactions_to_optimize) = obj_coef(2:5, 1);
                [solution] = calc_metabolic_metrics(constrained_model, [], ...
                    [], fva_grate, 'max', reactions_of_interest, [], 'fva');
                tissue_flux_values(match, :) = solution.maxflux;
            end
            
            [pearson_corr, pvalue] = corr(tissue_flux_values, ...
                tissue_matched_proteomics);
            reaction_names = reactions_of_interest(:, 3);
            
            tissue_structure = struct('Name', 'CCLE');
            fields = {...
                'Tissue', 'HistoneMark'; 'Reaction'; ...
                'PearsonR'; 'Pvalue'; 'FluxOutput'; 'ProteomicsOutput'
            };
        
            values = {...
                tissue, marks; reaction_names; ...
                pearson_corr; pvalue; tissue_flux_values; tissue_matched_proteomics
            };
        
            for i = 1:length(fields)
                tissue_structure.(fields{i}) = values{i};
            end
        end
end


 