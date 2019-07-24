%% Author: Marc Di Meo
clear;
initCobraToolbox(false);
changeCobraSolver('gurobi');

%% Load in files
load supplementary_software_code ...
    celllinenames_ccle1 ... % CCLE cell line names
    ccleids_met ... % Corresponding gene symbols
    ccle_expression_metz % Gene expression data (z-transformed)

path = './../new_var/';
vars = {...
    [path 'h3_ccle_names.mat'],... % CCLE cell line names for H3 proteomics data
    [path 'h3_marks.mat'],... % Correspond H3 histone marker IDs
    [path 'h3_media.mat'],... % Growth medium used in each experiment
    [path 'h3_relval.mat']... % Relative H3 proteomics data
    };

for kk = 1:numel(vars)
    load(vars{kk})
end

h3_relval = knnimpute(h3_relval);

%% Find common cell lines

common_celllines = intersect(h3_ccle_names, celllinenames_ccle1);
n = length(common_celllines);

%% Create cells for string arrays
h3_ccle_names_python = cellstr(h3_ccle_names);
h3_marks_python = cellstr(h3_marks);
common_celllines = cellstr(common_celllines);
save("correlation_value.mat")

%% Graphing correlation between
%plot(ccle_expression_metz, h3_relval)
