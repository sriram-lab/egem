initCobraToolbox
changeCobraSolver('gurobi');

addpath('/home/scampit/Desktop/eGEM/matlab/scripts/metabolic_sensitivity')
load ./../../metabolic_models/eGEM.mat
load ./../../vars/ccle_geneExpression_vars.mat
load ./../../vars/CCLE_epsilon.mat
reactions_of_interest = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
proteomics_path = './../../../data/CCLE/Global_Chromatin_Profiling/';
proteomics_array = readcell(strcat(proteomics_path, 'GCP_proteomics_remapped.csv'));

cell_names = proteomics_array(2:end, 1);
tissues = proteomics_array(2:end, 2);
marks = proteomics_array(1, 3:44);
medium = proteomics_array(2:end, 45);

proteomics = proteomics_array(2:end, 3:44);
missing_values = cellfun(@ismissing, proteomics);
proteomics(missing_values) = {NaN};
proteomics = cell2mat(proteomics);
proteomics = knnimpute(proteomics);

model = eGEM;
fva_grate = 100;
unique_tissues = unique(tissues);

for tiss = 1:length(unique_tissues)
    BIOMASS_OBJ_POS = find(ismember(eGEM.rxns, 'biomass_objective')); 

    tissue = unique_tissues(tiss);
    disp(tissue)
    tissue_positions = find(ismember(string(tissues), string(tissue)));
    tissue_proteomics = proteomics(tissue_positions, :);
    tissue_cellNames = cell_names(tissue_positions, :);
    tissue_medium = medium(tissue_positions, :);

    proteomics_CL_match_positions = find(ismember(string(tissue_cellNames), ...
        string(celllinenames_ccle1)));

    tissue_matched_proteomics = tissue_proteomics(proteomics_CL_match_positions, :);
    tissue_matched_cellNames = tissue_cellNames(proteomics_CL_match_positions,1);
    geneExp_CL_match_positions = find(ismember(string(celllinenames_ccle1),...
        string(tissue_matched_cellNames)));
    geneExp_CL_match_cellNames = celllinenames_ccle1(geneExp_CL_match_positions);
    [diffExp_genes] = find_diffexp_genes(eGEM, geneExp_CL_match_cellNames);

    for match = 1:length(proteomics_CL_match_positions)
        obj_coef = CCLE_epsilon.RPMI(2:5, 1);
        
        ON_fieldname = string(strcat('ON_', ...
            tissue_matched_cellNames(proteomics_CL_match_positions(match))));
        OFF_fieldname = string(strcat('OFF_', ...
            tissue_matched_cellNames(proteomics_CL_match_positions(match))));
        
        constrained_model = medium_LB_constraints(eGEM, tissue_medium(match));
        constrained_model.c(BIOMASS_OBJ_POS) = 1;
        
        reaction_name = char(reactions_of_interest(:, 1));
        reactions_to_optimize = [find(ismember(constrained_model.rxns,...
            reaction_name))];
        
        constrained_model.c(reactions_to_optimize) = obj_coef;
        
        kappa = 10;
        rho = 1E-3;
        epsilon = 1E-2;
        mode = 0; 
        epsilon2 = 1E-2;
        minfluxflag = true;
        [cellLine_model, solution] =  constrain_flux_regulation...
            (constrained_model, diffExp_genes.(ON_fieldname), ...
            diffExp_genes.(OFF_fieldname), ...
            [], [], [], mode, [], []);
        
        %[solution] = calc_metabolic_metrics(constrained_model, [], ...
        %    [], fva_grate, 'max', reactions_of_interest, [], 'fva');
        
        tissue_flux_values(match, :) = solution.flux(reactions_to_optimize);
    end
end

% [pearson_corr, pvalue] = corr(tissue_flux_values, ...
%     tissue_matched_proteomics);
% reaction_names = reactions_of_interest(:, 3);
% 
%     tissue_structure = struct('Name', 'CCLE');
%     fields = {...
%         'Tissue', 'HistoneMark'; 'Reaction'; ...
%         'PearsonR'; 'Pvalue'; 'FluxOutput'; 'ProteomicsOutput'
%     };
% 
%     values = {...
%         tissue, marks; reaction_names; ...
%         pearson_corr; pvalue; tissue_flux_values; tissue_matched_proteomics
%     };
% 
%     for i = 1:length(fields)
%         tissue_structure.(fields{i}) = values{i};
%     end
% end



 