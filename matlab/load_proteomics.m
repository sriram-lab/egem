function [cell_names, medium, marks, proteomics] = load_proteomics(proteomics_dataset_flag)
switch proteomics_dataset_flag
        case 'CCLE'
                load ccle_geneExpression_vars.mat
                load CCLE_Proteomics
                
        case 'LeRoy'
                load leroy.mat
                cell_names = leroy_cellline;
                marks = leroy_mark;
                medium = leroy_medium;
                proteomics = knnimpute(leroy_val);
                
end
end