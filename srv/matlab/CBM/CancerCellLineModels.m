%% Compute Cancer Cell Line COBRA models
% *Author:* Scott Campit
%% Summary
% This module computes metabolic fluxes for several histone post-translational 
% modification related demand reaction from a constrained eGEMM model using the 
% <https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results CCLE RNASeq> data 
% from <https://www.nature.com/articles/nature11003 Barretina et al., 2012>. 
% 
% *UPDATES:*
%% 
% * Included gamma score that takes into account the ratio between up-regulated 
% and down-regulated histone writers and erasers.
%% 1. Initialize COBRA
% Below are several preprocessig steps, including:
%% 
% * Setting biomass objective function
% * Lowering methionine substrate uptake rate 
% * Setting up params for iMAT (default settings used with pFBA)
% * Assigning the activity through demand reactions to be 0.0625. based on a 
% previous analysis with nutrient supplementation and deprivation.

% Initialize metabolic modeling components
clear all;
initCobraToolbox; changeCobraSolver('gurobi', 'all');
%% 2. Load metabolic models
% This code block loads specific versions of the eGEMM. There are currently 
% two versions of the eGEMM:
%% 
% # |07132020.mat:| contains all reactions from the metabolic model
% # |09072020_write_only.mat|: contains only histone writer reactions

load ~/Data/Reconstructions/eGEM/07132020.mat; % Contains all reactions, including demethylation and deacetylation reactions
eGEM_all = eGEM;
load ~/Data/Reconstructions/eGEM/09072020_write_only.mat
eGEM_write = eGEM;
models = {eGEM_all, eGEM_write};
%% 3. Add constraints to metabolic models
% The following code block sets the following constraints to the metabolic model:
%% 
% # Set the objective function to maximize biomass
% # Set the lower bound for methionine uptake
% # Let's turn off some specific reactions that cause trouble with metabolic 
% models.

for i = 1:length(models)
    mdl = models{i};
    
    % 1. Set up objective function
    mdl.c = zeros(size(mdl.c));
    mdl.c(string(mdl.rxns) == 'biomass_objective') = 1;
    
    % 2. Set methionine lower bound to a really small number
    [~, pos] = ismember({'EX_met_L_e'}, mdl.rxns);
    mdl.lb(pos) = -1;
    
    % 3. Turn off specific reactions
    posToTurnOff           = find(ismember(mdl.rxns, {'ALR', 'MGSA', 'MGSA2'}));
    mdl.ub(posToTurnOff) = 0; mdl.lb(posToTurnOff) = 0;
    
    % 4. Save the constrained model
    models{i} = mdl;
end
%% 4. Assign hyperparameters for the optimization procedure
% Assign default hyperparameter values for the iMAT algorithm.

% Set up hyperparameters for the linear programming section
hyperparams.eps = [];  hyperparams.isgenes = true;
hyperparams.kap = [];  hyperparams.eps2 = [];
hyperparams.rho = [];  hyperparams.pfba = false;
hyperparams.kap2 = []; hyperparams.gamma =[];
hyperparams.hscore = false;
%% 5. Set the acetylation and methylation objective coefficients to be a specific value
% Set the activity coefficient for demand reactions to be a fixed value, based 
% on the medium perturbation experiments I ran with this current iteration of 
% the eGEMM.

% Get list of demand reactions from eGEMM
pos = [3773:3777, 3780:3799];
dm_reactions = string(eGEM.rxns(pos))';

% Activity coefficient for demand reactions
epsilon_methylation = 0.069;
epsilon_acetylation = 0.069;

% Create an epsilon vector
for i = 1:length(dm_reactions)
    if contains(lower(dm_reactions(i)), "ac")
        epsilon(i, 1) = epsilon_acetylation;
    else
        epsilon(i, 1) = epsilon_methylation;
    end
end
%% 6. SANITY CHECK 1: Compute base metabolic model fluxes and growth rates for comparison
% Compute the solution for an unconstrained metabolic model. This will be the 
% control model.

for i = 1:length(dm_reactions)
    [~, rxnPos] = ismember(dm_reactions(i), eGEM.rxns);
    
    % Add activity to demand reactions
    tmp = eGEM;
    tmp.c(rxnPos) = epsilon(i);
    
    [~, soln] = CFR(tmp, hyperparams, {}, {});
    
    grate(1, i) = soln.x(3743);

    tmp2 = tmp;
    tmp2.lb(3743) = soln.x(3743) * 0.99;
    tmp2.c(3743)  = 0;                                                
    
    [~, soln] = CFR(tmp, hyperparams, {}, {});
    flux(1, i)  = soln.x(rxnPos);
end
%% 7. Preprocess Cancer Cell Line Encyclopedia data
% Load data from the Cancer Cell Line Encyclopedia
% This module uses RNASeq data from the CCLE and performs some basic preprocessing. 
% The data was normalized to TPM. There are 864 cancer cell lines by 1020 metabolic 
% genes that map onto RECON1 using Entrez identifiers.

load ~/Data/RNASeq/CCLE/CCLE_RNASeq.mat
load ~/Data/Proteomics/CCLE/CCLE_Proteomics.mat
% Get cell line intersections between the RNASeq dataset and the Global Chromatin Profiles

[~, ia, ib] = intersect(cell_lines, string(cell_names));
cell_lines = cell_lines(ia);
cell_names = cell_names(ib);

log2fc = log2fc(:, ia);
pvalue = pvalue(:, ia);

proteomics = proteomics(ib, :);
medium = medium(ib);
cultures = cultures(ib);
tissues = tissues(ib);
%% 8. Preprocess histone marker data
% Ensure that histone readers and writers are being used in the metabolic model
% First, let's get gene identifiers that are associated with histone writers 
% and erasers.

% Load up curated histone reaction map
histoneReactionFile = '~/Data/Reconstructions/eGEM/Epigenome-Scale Metabolic Reconstruction Map.xlsx';
histoneReactionMap = readcell(histoneReactionFile, 'Sheet', 'Genes');
histoneReactionMap(ismissing(string(histoneReactionMap(:, 8))), :) = [];

% Edit Entrez 
histone_entrez = cell2mat(histoneReactionMap(2:end, 8));
histone_entrez = string(histone_entrez);
% SANITY CHECK 2: Intersecting Entrez IDs between metabolic model and dataset
% This code block just checks and sees if there are any intersecting genes with 
% the metabolic model and the RNASeq data.

tmp1 = [];
tmp2 = [];
for j = 1:size(log2fc, 2)
        
    % Get differentially expressed genes for each cancer cell line.
    DE{j}.name    = cell_names(j);
    DE{j}.medium  = medium(j);
    DE{j}.culture = cultures(j);
    DE{j}.tissue  = tissues(j);
    DE{j}.upreg   = entrez_ids(log2fc(:, j) > 0 & pvalue(:, j) <= 0.05);
    DE{j}.downreg = entrez_ids(log2fc(:, j) < 0 & pvalue(:, j) <= 0.05);
    
    gene_list = [DE{j}.upreg; DE{j}.downreg];
    [~, ia, ib] = intersect(string(gene_list), string(histone_entrez));
    
    tmp1 = [tmp1; ia];
    tmp2 = [tmp2; ib];
    
end
%% 
% Since tmp1 and tmp2 are not empty, we are getting some histone markers.
%% 9. Compute CCLE metabolic fluxes from the eGEMM
% The sequence of steps for this code block is:
%% 
% # Iterate through reconstruction with all reactions and writers only
% # 

filepaths = [ ...
    "~/Data/CBM/eGEM/_Fluxes/CCLE_fluxes_all.mat", ...
    "~/Data/CBM/eGEM/_Fluxes/CCLE_fluxes_all_gamma.mat", ...
    "~/Data/CBM/eGEM/_Fluxes/CCLE_fluxes_writers.mat", ...
    "~/Data/CBM/eGEM/_Fluxes/CCLE_fluxes_writers_gamma.mat", ...
];

switches = ["No gamma", "Gamma"];

% Switch between reconstruction with all reactions and writers only
for i = 1:length(models)
    mdl = models{i};
    mdl.genes = mdl.geneEntrezID;
    entrez_ids = string(entrez_ids);
    
    flux = cell(length(dm_reactions), 1);
    grate = cell(length(dm_reactions), 1);
    
    % Run through all histone reactions
    for j = 1:length(dm_reactions)
        [~, rxnPos] = ismember(dm_reactions(j), mdl.rxns);
        tmp = mdl;
        tmp.c(rxnPos) = epsilon(j);
        
        % Run through all cancer cell lines
        for k = 1:size(log2fc, 2)
    
            % Get differentially expressed genes for each cancer cell line.
            DE{k}.name    = cell_names(k);
            DE{k}.medium  = medium(k);
            DE{k}.culture = cultures(k);
            DE{k}.tissue  = tissues(k);
            DE{k}.upreg   = entrez_ids(log2fc(:, k) > 0 & pvalue(:, k) <= 0.05);
            DE{k}.downreg = entrez_ids(log2fc(:, k) < 0 & pvalue(:, k) <= 0.05);
            
            % Constrain the metabolic model to be medium-specific
            mediumModel = addMediumConstraints(tmp, medium(k));
    
            % Switch between computing gamma (m == 2) and no gamma (m == 1)
            for m = 1:length(switches)
                
                if (i == 1 && m == 1) % All - gamma
                    hyperparams.gamma = [];
                    save(filepaths(1), 'flux', 'grate');
                elseif (i == 1 && m == 2) % All + gamma
                    hyperparams.gamma = compute_gamma(mdl, dm_reactions(j), histone_entrez, DE{k}.upreg, DE{k}.downreg);
                    save(filepaths(2), 'flux', 'grate');
                    
                elseif (i == 2 && m == 1) % Writers - gamma
                    hyperparams.gamma = [];
                    save(filepaths(3), 'flux', 'grate');
                    
                elseif (i == 2 && m == 2) % Writers + gamma
                    hyperparams.gamma = compute_gamma(mdl, dm_reactions(j), histone_entrez, DE{k}.upreg, DE{k}.downreg);
                    save(filepaths(4), 'flux', 'grate');
                end
                
                try
                    hyperparams.hscore = false;
                    % Compute flux using iMAT 
                    [mdl, soln] = CFR(mediumModel, hyperparams, ...
                                      DE{k}.upreg, DE{k}.downreg);
                                                   
                    % Get the growth rate
                    grates(k, 1) = soln.x(3743);
                    
                    % Ensure there's growth
                    tmp2 = mediumModel;
                    tmp2.lb(3743) = soln.x(3743) * 0.99;
                    tmp2.c(3743)  = 0;                                                
                    
                    if m == 2
                        hyperparams.hscore = true;
                    end
                    
                    % Compute metabolic flux when optimizing histone markers
                    [~, soln] = CFR(tmp2, hyperparams, DE{k}.upreg, DE{k}.downreg);
            
                    % Store data here
                    fluxes(:, k) = soln.x(1:length(eGEM.rxns));
                catch ME
                    grates(k, 1) = NaN;
                    fluxes(:, k) = NaN(length(eGEM.rxns), 1);
                end
            end
        end
        flux{j, 1} = fluxes;
        grate{j, 1} = grates;
        save(fileName, 'flux', 'grate', 'j', 'i', '-append');
    end 
end