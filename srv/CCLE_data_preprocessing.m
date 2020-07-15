%% CCLE Data Preprocessing
%% Summary
% This module preprocesses datasets for the cancer cell line encyclopedia.

clear all;
metabolomics = readcell('~/Data/Expression/Metabolomics/CCLE/CCLE_metabolomics.csv')
% Clean up the metabolomics dataset
% Next, we need to get the metabolomics dataset in a format that is appropriate 
% for creating maching learning models.
% 
% First, let's clean up some strings in the dataset by separating the tissue 
% with the cancer cell line.

% Save column names
colNames = metabolomics(1, 3:end);
colNames = ['Tissue', colNames];
metabolomics = metabolomics(2:end, [1, 3:end]);

% Separate by first '_'
[cellLine, remain] = strtok(metabolomics(:, 1)', '_');
tmp = strrep(remain, '_', ' '); cellLine = string(cellLine');
tissue = strtrim(tmp); tissue = string(tissue');

metabolomics(:, 1) = []; metabolomics = cell2mat(metabolomics);
% Find duplicates
% I took the first value for duplicate values.

[c, ia, idx] = unique(cellLine, 'stable');
metabolomics = metabolomics(ia, :);
tissue = tissue(ia); cellLine = cellLine(ia);

% Create table
% Now we'll save the data as a table for future reference.

% Construct table
metabolomics = array2table([tissue, metabolomics], ...
    "RowNames", string(cellLine)', ...
    "VariableNames", colNames);

% Save as file
save('ccle_metabolomics.mat', 'metabolomics');
% Clean up the GCP dataset
% Find duplicates
% I took the first value for duplicate values.

% Map cell lines to proteomics data
[c, ia, idx] = unique(cell_names, 'stable');
cell_names = string(cell_names); cell_names = cell_names(ia);
proteomics = proteomics(ia, :);
% Create table
% Now we'll save the data as a table for future reference.

% Construct table
gcp = array2table(proteomics, ...
    "RowNames", cell_names, ...
    "VariableNames", marks);

% Save as file
save('CCLE_Proteomics.mat', '-append', 'gcp');
% Clean up DNA Microarray dataset

load ccle_geneExpression_vars.mat
% Find duplicates
% I took the first value for duplicate values.

% Get unique gene symbols
[c, ia, idx] = unique(ccleids_met, 'stable');
ccleids_met = string(ccleids_met); ccleids_met = ccleids_met(ia);
ccle_expression_metz = ccle_expression_metz(ia, :);

% Get unique cell line names
[c, ia, idx] = unique(celllinenames_ccle1, 'stable');
celllinenames_ccle1 = string(celllinenames_ccle1); 
celllinenames_ccle1 = celllinenames_ccle1(ia);
ccle_expression_metz = ccle_expression_metz(:, ia);
% Create table
% Now we'll save the data as a table for future reference.

% Construct table
microarray = array2table(ccle_expression_metz, ...
    "RowNames", ccleids_met, ...
    "VariableNames", celllinenames_ccle1);

% Save as file
save('ccle_geneExpression_vars.mat', '-append', 'microarray');
% Clean up RNASeq dataset
% Next, we need to get the metabolomics dataset in a format that is appropriate 
% for creating maching learning models.

clear all;
metabolomics = readcell('~/Data/Expression/RNASeq/CCLE/CCLE_metabolomics.csv')
%% 
% First, let's clean up some strings in the dataset by separating the tissue 
% with the cancer cell line.

% Save column names
colNames = metabolomics(1, 3:end);
colNames = ['Tissue', colNames];
metabolomics = metabolomics(2:end, [1, 3:end]);

% Separate by first '_'
[cellLine, remain] = strtok(metabolomics(:, 1)', '_');
tmp = strrep(remain, '_', ' '); cellLine = string(cellLine');
tissue = strtrim(tmp); tissue = string(tissue');

metabolomics(:, 1) = []; metabolomics = cell2mat(metabolomics);
% Find duplicates
% I took the first value for duplicate values.

[c, ia, idx] = unique(cellLine, 'stable');
metabolomics = metabolomics(ia, :);
tissue = tissue(ia); cellLine = cellLine(ia);

% Create table
% Now we'll save the data as a table for future reference.

% Construct table
metabolomics = array2table([tissue, metabolomics], ...
    "RowNames", string(cellLine)', ...
    "VariableNames", colNames);

% Save as file
save('ccle_metabolomics.mat', 'metabolomics');