%% Author: Marc Di Meo

clear;

initCobraToolbox(false);
changeCobraSolver('gurobi');


%% Load in files
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz % contains CCLE cellline names for gene exp, enzymes encoding specific metabolites, and gene expression data (z-transformed)



path = './../new_var/';
vars = {...
    [path 'h3_ccle_names.mat'], [path 'h3_marks.mat'],...
    [path 'h3_media.mat'], [path 'h3_relval.mat']...
    }; % contains CCLE cellline names for H3 proteomics, corresponding marker ids, growth media, relative H3 proteomics
for kk = 1:numel(vars)
    load(vars{kk})
end

h3_relval = knnimpute(h3_relval); 

%% Find common cell lines

common_celllines = intersect(h3_ccle_names, celllinenames_ccle1);
n = length(common_celllines); %used for Pearson Correlation

%% Create cells for string arrays

h3_ccle_names_python = cellstr(h3_ccle_names);
h3_marks_python = cellstr(h3_marks);
common_celllines = cellstr(common_celllines);

save('Recon1_genes.mat', RECON1.genes)
save("correlation_value.mat")

%% Fooling around




%% Graphing correlation between
%plot(ccle_expression_metz, h3_relval)



