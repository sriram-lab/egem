function [cell_names, medium, marks, proteomics] = load_proteomics(proteomics_dataset_flag)
    switch proteomics_dataset_flag
        case 'CCLE'
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
            proteomics = knnimpute(proteomics);

    end
end