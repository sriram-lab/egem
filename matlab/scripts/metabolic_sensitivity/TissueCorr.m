function [tissue_structure] = TissueCorr(tissue_flux_values, ...
<<<<<<< HEAD
    tissue_matched_proteomics)

    reaction_names = {'Histone 1Me', 'Histone 2Me', ...
        'Histone 3Me', 'Histone Ac'};
=======
    tissue_matched_proteomics, tissue, marks)

    reaction_names = {'Histone Ac', 'Histone 1Me', 'Histone 2Me', ...
        'Histone 3Me'};
>>>>>>> 4ad2bfc7bcf688e7fc426cfe727c05970b7421f5
    
    [pearson_corr, pvalue] = corr(tissue_flux_values, ...
        tissue_matched_proteomics);

    tissue_structure = struct('Name', 'CCLE');
    fields = {...
<<<<<<< HEAD
        'Tissue', 'HistoneMark'; 'Reaction'; ...
        'PearsonR'; 'Pvalue'; 'FluxOutput'; 'ProteomicsOutput'
    };
    values = {...
        tissue, marks; reaction_names; ...
=======
        'Tissue'; 'HistoneMark'; 'Reaction'; ...
        'PearsonR'; 'Pvalue'; 'FluxOutput'; 'ProteomicsOutput'
    };
    values = {...
        tissue; marks; reaction_names; ...
>>>>>>> 4ad2bfc7bcf688e7fc426cfe727c05970b7421f5
        pearson_corr; pvalue; tissue_flux_values; tissue_matched_proteomics
    };

    for i = 1:length(fields)
        tissue_structure.(fields{i}) = values{i};
    end
end