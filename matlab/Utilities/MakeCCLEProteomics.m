% Make histone proteomics variable

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

save('CCLE_Proteomics.mat', 'proteomics', 'cell_names', 'tissues', 'marks', 'medium');
