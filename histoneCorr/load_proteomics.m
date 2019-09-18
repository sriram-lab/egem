function [cell_names, medium, marks, proteomics] = load_proteomics(proteomics_dataset_flag)
    switch proteomics_dataset_flag
        case 'CCLE'
            load ccle_geneExpression_vars.mat
            load CCLE_Proteomics

        case 'LeRoy'
            vars = {...
                ['leroy_cellline.mat'],... % CCLE cellline names for H3 proteomics, 
                ['leroy_mark.mat'],... % Marker IDs
                ['leroy_val.mat'],...% Average values
                ['leroy_medium.mat']
                }; 

            for kk = 1:numel(vars) 
                load(vars{kk})
            end

            cell_names = cell(:,1);
            marks = leroy_mark;

            proteomics = leroy_val;
            proteomics = proteomics';
            proteomics = knnimpute(proteomics);

    end
end