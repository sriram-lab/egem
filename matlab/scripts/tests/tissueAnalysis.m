%initCobraToolbox
%changeCobraSolver('gurobi');
function tissueAnalysis()
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
unique_tissues = unique(tissues);
BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
model.c(BIOMASS_OBJ_POS) = 1;

for tiss = 1:length(unique_tissues)
    oneTissue = unique_tissues(tiss);
    disp(oneTissue)
    
    tissuePositions = find(ismember(string(tissues), string(oneTissue)));
    tissueCellNames = cell_names(tissuePositions, :);
    
    tissueToGeneMatch = find(ismember(string(celllinenames_ccle1), ...
        string(tissueCellNames)));
    geneExpMatch = celllinenames_ccle1(tissueToGeneMatch);

    tissueProteomicsValues = proteomics(tissuePositions, :);
    tissueMedium = medium(tissuePositions, :);
    
    [diffExp_genes] = find_diffexp_genes(model, geneExpMatch);

    for match = 1:length(tissueCellNames)
        obj_coef = [1E-3, 1E-3, 1E-3, 1E-3];

        ON_fieldname = string(strcat('ON_', tissueCellNames(match)));
        OFF_fieldname = string(strcat('OFF_', tissueCellNames(match)));

        constrained_model = medium_LB_constraints(model, tissueMedium(match));

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
            
            tissue_flux_values(match, rxn) = solution.flux(optimized_rxn);
            tissue_grates(match, rxn) = solution.flux(BIOMASS_OBJ_POS);
            constrained_model.c(reactions_to_optimize) = 0;
            
            tissue = string(oneTissue);
            tissue = strrep(tissue, ' ', '');
            fluxVarName = string(strcat(tissue, '_flux_values(match, rxn) = solution.flux(optimized_rxn);'));
            grateVarName = string(strcat(tissue, '_grates(match, rxn) = solution.x(BIOMASS_OBJ_POS);'));
            eval(fluxVarName);
            eval(grateVarName);
        end
    end
    
    corrVarName = string(strcat(tissue, '_corr = ', 'TissueCorr(', ...
        tissue, '_flux_values, tissueProteomicsValues, tissue, marks)'));
    eval(corrVarName)
    
    plotVarName = string(strcat('plot_heatmap(', tissue, ...
        '_corr', ', [], "correlation", [], [], [])'));
    eval(plotVarName)
end

 