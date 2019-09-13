function [solution] = mediaAnalysis(exp)
    %addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\metabolic_sensitivity')
    %addpath('C:\Users\scampit\Desktop\egem\matlab\scripts\visualizations')
    
    addpath('/home/scampit/Desktop/egem/matlab/scripts/nutrient_sensitivity')
    addpath('/home/scampit/Desktop/egem/matlab/scripts/visualizations')
    
    load ./../../metabolic_models/eGEM_mm.mat
    load ./../../vars/ccle_geneExpression_vars.mat
    load ./../../vars/CCLE_Proteomics
    load allVars.mat
    
    kappa = 1;
    rho = 1;
    epsilon = 1E-3;
    mode = 0; 
    epsilon2 = 1E-3;
    minfluxflag = true;

    reactions_of_interest = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
    [~, mediaList] = xlsfinfo('./../../../data/Medium_Component_Maps/final_medium2.xlsx');
    model = eGEM_mm;
    BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); 
    model.c(BIOMASS_OBJ_POS) = 1;
    unique_medium = unique(medium);
    mediaList = intersect(string(unique_medium), string(mediaList));

    for med = 1:length(mediaList)
        disp(mediaList(med))
        medium = strtrim(medium);
        mediaPositions = find(ismember(string(medium), string(mediaList(med))));
        mediaCellNames = cell_names(mediaPositions, :);

        mediaToGeneMatch = find(ismember(string(celllinenames_ccle1), ...
            string(mediaCellNames)));
        geneExpMatch = celllinenames_ccle1(mediaToGeneMatch);

        matched = intersect(string(mediaCellNames), string(geneExpMatch));
        matchedMedia = find(ismember(string(cell_names), string(matched)));
        matchedGenes = find(ismember(string(celllinenames_ccle1), string(matched)));
        mediaCellNames = cell_names(matchedMedia, :);
        geneExpMatch = celllinenames_ccle1(matchedGenes);

        mediaProteomicsValues = proteomics(matchedMedia, :);

        [diffExp_genes] = find_diffexp_genes(model, geneExpMatch);

        for match = 1:length(mediaCellNames)

            ON_fieldname = string(strcat('ON_', geneExpMatch(match)));
            OFF_fieldname = string(strcat('OFF_', geneExpMatch(match)));
            
            str = strcat(string("ObjCoef = epsilon2_"), string(lower(mediaList(med))), ...
                string("(:, 1)"));
            eval(str);

            constrained_model = medium_LB_constraints(model, mediaList(med));

            reaction_name = char(reactions_of_interest(:, 1));
            reactions_to_optimize = [find(ismember(constrained_model.rxns,...
                reaction_name))];
            
            final_medium = string(mediaList(med));
            final_medium = strrep(final_medium, '-', '');
            final_medium = strrep(final_medium, ' ', '');

            switch exp
                case {'SRA', 'NoComp'}
                    for rxn = 1:length(reactions_to_optimize)
                        optimized_rxn = reactions_to_optimize(rxn);
                        constrained_model.c(reactions_to_optimize) = ObjCoef(rxn);

                        [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);
                        
                        fluxVarName = string(strcat(final_medium, exp, 'flux(match, rxn) = solution.flux(optimized_rxn);'));
                        grateVarName = string(strcat(final_medium, exp, 'grates(match, rxn) = solution.x(BIOMASS_OBJ_POS);'));
                        eval(fluxVarName);
                        eval(grateVarName);
                        constrained_model.c(reactions_to_optimize) = 0;
                        
                    end
                    
                case 'Comp'
                    reaction_positions = find(ismember(model.rxns, reactions_of_interest));
                    constrained_model.c(reactions_to_optimize) = ObjCoef;

                    [~, solution] =  constrain_flux_regulation...
                           (constrained_model, diffExp_genes.(ON_fieldname), ...
                            diffExp_genes.(OFF_fieldname), ...
                            kappa, rho, epsilon, mode, [], minfluxflag);

                    fluxVarName = string(strcat(final_medium, exp, 'flux(match, :) = solution.flux(reaction_positions);'));
                    grateVarName = string(strcat(final_medium, exp, 'grates(match, 1) = solution.x(BIOMASS_OBJ_POS);'));
                    eval(fluxVarName);
                    eval(grateVarName);
                    
                case 'FVA'
                    reaction_positions = find(ismember(model.rxns, reactions_of_interest));
                    constrained_model.c(reaction_positions) = ObjCoef;
                    [~, maxFlux] = fluxVariability(constrained_model, 100, ...
                            'max', reactions_of_interest);
                    fluxVarName = string(strcat(final_medium, exp, 'flux(match, :) = maxFlux;'));
                    eval(fluxVarName);                    
    
            end
        end

        corrVarName = string(strcat(final_medium, exp, 'corr = ', 'MediumCorr(', ...
            final_medium, exp, 'flux, mediaProteomicsValues, mediaList(med), marks)'));
        eval(corrVarName)

        plotVarName = string(strcat('plotHistoneCorrelation(', final_medium, ...
            exp, 'corr, "correlation", exp, "medium_corr")'));
        eval(plotVarName)
        
        BPVarName = string(strcat("makeBoxPlot(", final_medium, exp, ...
            "flux, mediaProteomicsValues, final_medium, 'medium_bp', exp)"));
        eval(BPVarName)
        
        histVarName = string(strcat("makeHist(", final_medium, exp, ...
            "flux, string(mediaList(med)), 'medium_hist', exp)"));
        eval(histVarName)
        
        fluxFileName = './../../vars/histoneFluxValues.mat';
        if exist(fluxFileName, 'file')
            saveStr1 = strcat("save(fluxFileName, '", ...
                final_medium, exp, "flux', '-append')");
            eval(string(saveStr1)); 
        else
            saveStr1 = strcat("save(fluxFileName, '", ...
                final_medium, exp, "flux')");
            eval(string(saveStr1));
        end
        
        corrFileName = './../../vars/histoneCorrValues.mat';
        if exist(corrFileName, 'file')
            saveStr2 = strcat("save(corrFileName, '", ...
                final_medium, exp, "corr', '-append')");
            eval(string(saveStr2));
        else
            saveStr2 = strcat("save(corrFileName, '", ...
                final_medium, exp, "corr')");
            eval(string(saveStr2));
        end
    end
end

 