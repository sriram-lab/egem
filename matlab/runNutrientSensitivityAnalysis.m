initCobraToolbox;
changeCobraSolver('gurobi');

%% Optimization 1A: Run Single reaction activity (SRA) for histone reactions
load ./../../vars/CCLE_Proteomics medium
load ./../../metabolic_models/eGEM_mm
model = eGEM_mm;
[~, mediaList] = xlsfinfo('./../../../data/Medium_Component_Maps/final_medium2.xlsx');
unique_medium = unique(mediaList);
mediaList = intersect(string(unique_medium), string(mediaList));

epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for medium = 1:length(mediaList)
    disp(mediaList(medium))
    for objCoef = 1:length(epsilon2)
        sraStr =  strcat("[sra", string(objCoef), '_', mediaList(medium),...
            "] = nutrient_sensitivity(model, mediaList(medium), ", ...
            "epsilon2(objCoef), 'sra')");
        eval(sraStr);
    end
end

%% Calculate epsilon2 values to use for fba by dynamic range
for i=1:length(mediaList)
    str = strcat("epsilon2_",lower(mediaList(i)), " = ", ...
        "dynamic_range(sra1_", mediaList(i), ", ", ...
        "sra2_", mediaList(i), ", ", ...
        "sra3_", mediaList(i), ", ", ...
        "sra4_", mediaList(i), ", ", ...
        "sra5_", mediaList(i), ", ", ...
        "sra6_", mediaList(i), ", ", ...
        "sra7_", mediaList(i), ", ", ...
        "'dynamic');");
    eval(str);
end

%% Optimization procedure using no competition and competition
for medium = 1:length(mediaList)
    disp(mediaList(medium))
    
    noCompStr =  strcat("[fba_", lower(mediaList(medium)),"_noCompetition]", ...
        "= nutrient_sensitivity(model, mediaList(medium), epsilon2_", ...
        lower(mediaList(medium)), ", 'no_competition')");
    eval(noCompStr);
    
    noCompPlot = strcat("plotNutrientSensitivity(fba_", lower(mediaList(medium)),...
        "_noCompetition, 'no_competition', epsilon2, mediaList(medium))");
    eval(noCompPlot);
end

for medium = 1:length(mediaList)
    disp(mediaList(medium))

    CompStr =  strcat("[fba_", lower(mediaList(medium)),"_Competition]", ...
        "= nutrient_sensitivity(model, mediaList(medium), epsilon2_", ...
        lower(mediaList(medium)), ", 'competition')");
    eval(CompStr);
    
    CompPlot = strcat("plotNutrientSensitivity(fba_", lower(mediaList(medium)),...
        "_Competition, 'competition', epsilon2, mediaList(medium))");
    eval(CompPlot);
end

for medium = 1:length(mediaList)
    disp(mediaList(medium))

    fvaStr =  strcat("[fva_", lower(mediaList(medium)),"_fva]", ...
        "= nutrient_sensitivity(model, mediaList(medium), epsilon2_", ...
        lower(mediaList(medium)), ", 'fva')");
    eval(fvaStr);
    
    fvaPlot = strcat("plotNutrientSensitivity(fva_", lower(mediaList(medium)),...
        "_fva, 'fva', epsilon2, mediaList(medium))");
    eval(fvaPlot);
end