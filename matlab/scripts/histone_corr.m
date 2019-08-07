%% histone_corr 
% @author: Scott Campit
function [STRUCT] = histone_corr(model, reactions_of_interest,...
    eps_struct, mode, epsilon, rho, kappa, minfluxflag, exp,...
    proteomics_dataset, fva_grate)

switch proteomics_dataset
    case 'CCLE'
        % Load the relative H3 proteomics dataset from the CCLE
        path = './../data/';
        proteomics_array = readcell(strcat(path, 'GCP_proteomics_remapped.csv'));
       
        cell_names = proteomics_array(2:end, 1);
        tissues = proteomics_array(2:end, 2);
        marks = proteomics_array(1, 3:44);
        medium = proteomics_array(2:end, 45);
        proteomics = proteomics_array(2:end, 3:44);
        
    case 'LeRoy'
        path = './../new_var/';
        vars = {...
            [path1 'leroy_cellline.mat'],... % CCLE cellline names for H3 proteomics, 
            [path1 'leroy_mark.mat'],... % Marker IDs
            [path1 'leroy_val.mat'],...% Average values
            }; 
        
        for kk = 1:numel(vars) 
            load(vars{kk})
        end
        
        cell_names = cell(:,1);
        medium = cell(:,2);
        marks = leroy_mark;
        proteomics = leroy_val;
        proteomics = proteomics';
end
        
load ./../vars/supplementary_software_code...
    celllinenames_ccle1... % CCLE cellline names
    ccleids_met... % Gene symbols
    ccle_expression_metz % Z-transformed gene expression

BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
proteomics = knnimpute(proteomics);
proteomics_CL_match = find(ismember(cell_names, celllinenames_ccle1));
number_of_matched_proteomics = length(proteomics_CL_match);
proteomics = proteomics(proteomics_CL_match, :);
cell_names = cell_names(proteomics_CL_match,1);
geneExp_CL_match = find(ismember(celllinenames_ccle1, cell_names));

for match = 1:number_of_matched_proteomics
    disp(['Cell line: ', cell_names(match)])
    constrained_model = model;
    
    ongenes = unique(ccleids_met(ccle_expression_metz(:, geneExp_CL_match(match)) >= 2));
    offgenes = unique(ccleids_met(ccle_expression_metz(:, geneExp_CL_match(match)) <= -2));
    ongenes = intersect(ongenes, constrained_model.rxns);
    offgenes = intersect(offgenes, constrained_model.rxns);
    [ix, pos]  = ismember({'EX_met_L(e)'}, constrained_model.rxns);
    
    disp(medium(match))
    obj_coef = eps_struct.RPMI; % This is temporary. I will have to recalculate these values based on the new formulation.
    
    constrained_model.lb(pos) = -0.5;
    constrained_model = media(constrained_model, medium(match));
    constrained_model.c(BIOMASS_OBJ_POS) = 1;

    [~,~,onreactions,~] =  deleteModelGenes(constrained_model, ongenes);
    [~,~,offreactions,~] =  deleteModelGenes(constrained_model, offgenes);

    rxnname = char(reactions_of_interest(:, 1));
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
            rxnpos = [find(ismember(model_to_test.rxns, rxnname))];
            model_to_test.c(rxnpos) = obj_coef(:, 1);
            [flux, ~, ~] =  constrain_flux_regulation(model_to_test,...
                onreactions, offreactions, kappa, rho, epsilon, mode , [], ...
                minfluxflag);
            all_flux_values(match,:) = flux(rxnpos);
        case 'fva'
            model_to_test = constrained_model;
            rxnpos = [find(ismember(model_to_test.rxns, rxnname))];
            model_to_test.c(rxnpos) = obj_coef(:, 1);
            [~, ~, ~, ~, flux, ~] =...
                calc_metabolic_metrics(model_to_test, [], [], fva_grate,...
                'max', reactions_of_interest, [], 'fva');
            all_flux_values(match,:) = flux;
    end
end

[rho, pval] = corr(all_flux_values, proteomics);
rxns = reactions_of_interest(:, 3);

%% Save data in struct
STRUCT = struct('Name', proteomics_dataset);

fields = {...
        'HistoneMark'; 'Reaction'; ...
        'R'; 'Pvalue'; 'Flux'; 'Proteomics'
    };

values = {...
    marks; rxns; ...
    rho; pval; all_flux_values; proteomics
    };

for match=1:length(fields)
    STRUCT.(fields{match}) = values{match};
end

end