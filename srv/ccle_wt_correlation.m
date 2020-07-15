%% Correlation analysis from wild type CCLE models
%% Summary
% 
%% Load data and get intersection between transcriptomics and GCP

load ccle_geneExpression_vars.mat
load CCLE_Proteomics.mat

% Get intersection indicies
[~, ia, ib] = intersect(string(celllinenames_ccle1), string(cell_names));

% Extract intersection data
microarray_cellnames = celllinenames_ccle1(ia);
microarray_data = ccle_expression_metz(:, ia);
microarray_genenames = ccleids_met;

gcp_cellnames = cell_names(ib);
gcp_medium = medium(ib);
gcp_data = proteomics(ib, :);
gcp_culture = cultures(ib, :);
gcp_tissue = tissues(ib, :);
%% Load flux data from single reaction optimization and get cell lines that have a growth rate

load ccle_flux_sra.mat

% Turn table to array to make processing easier
colnames = histone_flux.Properties.VariableNames;
rowNames = histone_flux.Properties.RowNames;
histone_flux = table2array(histone_flux);

% Get cell lines with growth rate
hist_data = histone_flux(histone_flux(:, 3743) > 0, :)
cell_lines = rowNames(histone_flux(:, 3743) > 0);
filtered_gcp = rescale(gcp_data(histone_flux(:, 3743) > 0, :));
%% Get correlation from all cancer cell lines that have a growth rate 
% I selected specific histone markers to analyze that are known to respond to 
% nutrient perturbation. For the most part, there's no correlation

% Manually check specific correlation between histone markers and bulk
% histone levels
corr(hist_data(:, 3756), filtered_gcp(:, 2), 'Type', 'Pearson')  % H3K4me1
corr(hist_data(:, 3757), filtered_gcp(:, 3), 'Type', 'Pearson')  % H3K4me2
corr(hist_data(:, 3755), filtered_gcp(:, 4), 'Type', 'Pearson')  % H3K4ac1

corr(hist_data(:, 3756), filtered_gcp(:, 6), 'Type', 'Pearson')  % H3K9me1
corr(hist_data(:, 3757), filtered_gcp(:, 7), 'Type', 'Pearson')  % H3K9me2
corr(hist_data(:, 3758), filtered_gcp(:, 8), 'Type', 'Pearson')  % H3K9me3
corr(hist_data(:, 3755), filtered_gcp(:, 9), 'Type', 'Pearson')  % H3K9ac1

corr(hist_data(:, 3756), filtered_gcp(:, 24), 'Type', 'Pearson') % H3K27me1
corr(hist_data(:, 3757), filtered_gcp(:, 28), 'Type', 'Pearson') % H3K27me2
corr(hist_data(:, 3757), filtered_gcp(:, 31), 'Type', 'Pearson') % H3K27me3
corr(hist_data(:, 3755), filtered_gcp(:, 33), 'Type', 'Pearson') % H3K27ac1

corr(hist_data(:, 3756), filtered_gcp(:, 21), 'Type', 'Pearson') % H3K36me1
corr(hist_data(:, 3757), filtered_gcp(:, 22), 'Type', 'Pearson') % H3K36me2
corr(hist_data(:, 3758), filtered_gcp(:, 23), 'Type', 'Pearson') % H3K36me3

corr(hist_data(:, 3756), filtered_gcp(:, 41), 'Type', 'Pearson') % H3K79me1
corr(hist_data(:, 3757), filtered_gcp(:, 42), 'Type', 'Pearson') % H3K79me2
%% Get correlation from cancer cell lines that are from the TCGA-110 CCL Panel

% Get TCGA-110 CCL
load tcga_meta.mat
tcga_cellname = tcga_meta(:, 'CellLine');
[names, ia, ib] = intersect(cell_lines, string(table2array(tcga_cellname)));

% Manually check specific correlation between histone markers and bulk
% histone levels
corr(hist_data(ia, 3756), filtered_gcp(ia, 2), 'Type', 'Pearson')  % H3K4me1
corr(hist_data(ia, 3757), filtered_gcp(ia, 3), 'Type', 'Pearson')  % H3K4me2
corr(hist_data(ia, 3755), filtered_gcp(ia, 4), 'Type', 'Pearson')  % H3K4ac1

corr(hist_data(ia, 3756), filtered_gcp(ia, 6), 'Type', 'Pearson')  % H3K9me1
corr(hist_data(ia, 3757), filtered_gcp(ia, 7), 'Type', 'Pearson')  % H3K9me2
corr(hist_data(ia, 3758), filtered_gcp(ia, 8), 'Type', 'Pearson')  % H3K9me3
corr(hist_data(ia, 3755), filtered_gcp(ia, 9), 'Type', 'Pearson')  % H3K9ac1

corr(hist_data(ia, 3756), filtered_gcp(ia, 24), 'Type', 'Pearson') % H3K27me1
corr(hist_data(ia, 3757), filtered_gcp(ia, 28), 'Type', 'Pearson') % H3K27me2
corr(hist_data(ia, 3757), filtered_gcp(ia, 31), 'Type', 'Pearson') % H3K27me3
corr(hist_data(ia, 3755), filtered_gcp(ia, 33), 'Type', 'Pearson') % H3K27ac1

corr(hist_data(ia, 3756), filtered_gcp(ia, 21), 'Type', 'Pearson') % H3K36me1
corr(hist_data(ia, 3757), filtered_gcp(ia, 22), 'Type', 'Pearson') % H3K36me2
corr(hist_data(ia, 3758), filtered_gcp(ia, 23), 'Type', 'Pearson') % H3K36me3

corr(hist_data(ia, 3756), filtered_gcp(ia, 41), 'Type', 'Pearson') % H3K79me1
corr(hist_data(ia, 3757), filtered_gcp(ia, 42), 'Type', 'Pearson') % H3K79me2
%% 
%