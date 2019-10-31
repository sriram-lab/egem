function [solution] = tissueAnalysis(exp)
    
    addpath('/home/scampit/Desktop/egem/matlab/scripts/nutrient_sensitivity')
    addpath('/home/scampit/Desktop/egem/matlab/scripts/visualizations')
    
    load ./../../metabolic_models/eGEM_mm.mat
    load ./../../vars/ccle_geneExpression_vars.mat
    load ./../../vars/CCLE_Proteomics
    load allVars.mat

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

        matched = intersect(string(tissueCellNames), string(geneExpMatch));
        matchedTissue = find(ismember(string(cell_names), string(matched)));
        matchedGenes = find(ismember(string(celllinenames_ccle1), string(matched)));
        tissueCellNames = cell_names(matchedTissue, :);
        geneExpMatch = celllinenames_ccle1(matchedGenes);

        tissueProteomicsValues = proteomics(matchedTissue, :);
        tissueMedium = medium(matchedTissue, :);

        [diffExp_genes] = find_diffexp_genes(model, geneExpMatch);

        for match = 1:length(tissueCellNames)

            ON_fieldname = string(strcat('ON_', geneExpMatch(match)));
            OFF_fieldname = string(strcat('OFF_', geneExpMatch(match)));
            str = strcat(string("ObjCoef = epsilon2_"), string(lower(tissueMedium(match))), string("(:, 1)"));
            eval(str);            
            
            constrained_model = medium_LB_constraints(model, tissueMedium(match));

            reaction_name = char(reactions_of_interest(:, 1));
            reactions_to_optimize = [find(ismember(constrained_model.rxns,...
                reaction_name))];
            constrained_model.c(reactions_to_optimize) = 0;
            
            tissue = string(oneTissue);
            tissue = strrep(tissue, ' ', '');
            
            kappa = 1;
            rho = 1;
            epsilon = 1E-3;
            mode = 0; 
            epsilon2 = 1E-3;
            minfluxflag = true;
            
            switch exp
                case {'SRA', 'NoComp'}
                    for rxn = 1:length(reactions_to_optimize)
                        optimized_rxn = reactions_to_optimize(rxn);
                        constrained_model.c(reactions_to_optimize) = ObjCoef(rxn);

                        [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);
                        
                        fluxVarName = string(strcat(tissue, exp, 'flux(match, rxn) = solution.flux(optimized_rxn);'));
                        grateVarName = string(strcat(tissue, exp, 'grates(match, rxn) = solution.x(BIOMASS_OBJ_POS);'));
                        eval(fluxVarName);
                        eval(grateVarName);
                        constrained_model.c(reactions_to_optimize) = 0;         
                    end
                    
                case 'Comp'
                    constrained_model.c(reactions_to_optimize) = ObjCoef;

                    [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);

                    fluxVarName = string(strcat(tissue, exp, 'flux(match, :) = solution.flux(reactions_to_optimize);'));
                    grateVarName = string(strcat(tissue, exp, 'grates(match, 1) = solution.x(BIOMASS_OBJ_POS);'));
                    eval(fluxVarName);
                    eval(grateVarName);
                    
                case 'FVA'
                    model.c(reactions_to_optimize) = ObjCoef;
                    [~, maxFlux] = fluxVariability(constrained_model, 100, ...
                            'max', reactions_of_interest);
                    fluxVarName = string(strcat(tissue, exp, 'flux(match, :) = maxFlux;'));
                    eval(fluxVarName);                       
            end
        end
        
        % Make plots and save data
        corrVarName = string(strcat(tissue, exp, 'corr = ', 'TissueCorr(', ...
           tissue, exp, 'flux, tissueProteomicsValues, tissue, marks)'));
        eval(corrVarName)

        plotVarName = string(strcat('plotHistoneCorrelation(', tissue, exp, ...
           'corr, "correlation", exp, "tissue_corr")'));
        eval(plotVarName)
        
        BPVarName = string(strcat("makeBoxPlot(", tissue, exp, ...
            "flux, tissueProteomicsValues, string(oneTissue), 'tissue_bp', exp)"));
        eval(BPVarName)
        
        histVarName = string(strcat("makeHist(", tissue, exp, ...
            "flux, string(oneTissue), 'tissue_hist', exp)"));
        eval(histVarName)
        
        fluxFileName = './../../vars/histoneFluxValues.mat';
        if exist(fluxFileName, 'file')
            saveStr1 = strcat("save(fluxFileName, '", ...
                tissue, exp, "flux', '-append')");
            eval(string(saveStr1)); 
        else
            saveStr1 = strcat("save(fluxFileName, '", ...
                tissue, exp, "flux')");
            eval(string(saveStr1));
        end
        
        corrFileName = './../../vars/histoneCorrValues.mat';
        if exist(corrFileName, 'file')
            saveStr2 = strcat("save(corrFileName, '", ...
                tissue, exp, "corr', '-append')");
            eval(string(saveStr2));
        else
            saveStr2 = strcat("save(corrFileName, '", ...
                tissue, exp, "corr')");
            eval(string(saveStr2));
        end      
    end
end

 