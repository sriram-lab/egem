%% Compute Cancer Cell Line COBRA models
% *Author:* Scott Campit
%% Summary
% This module computes metabolic fluxes for several histone post-translational 
% modification related demand reaction from a constrained eGEMM model using the 
% <https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results CCLE RNASeq> data 
% from <https://www.nature.com/articles/nature11003 Barretina et al., 2012>. 
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

load /nfs/turbo/umms-csriram/scampit/Data/Reconstructions/eGEM/07132020.mat; % Contains all reactions, including demethylation and deacetylation reactions
%load /nfs/turbo/umms-csriram/scampit/Data/Reconstructions/eGEM/09072020_write_only.mat
%% 
% Set the objective function to maximize biomass

% Set up objective function
eGEM.c = zeros(size(eGEM.c));
eGEM.c(string(eGEM.rxns) == 'biomass_objective') = 1;
%% 
% Set the lower bound for methionine uptake

% Set methionine lower bound to a really small number
[~, pos] = ismember({'EX_met_L_e'}, eGEM.rxns);
eGEM.lb(pos) = -1;
%% 
% Let's turn off some specific reactions that cause trouble with metabolic models.

% Turn off specific reactions
posToTurnOff           = find(ismember(eGEM.rxns, {'ALR', 'MGSA', 'MGSA2'}));
eGEM.ub(posToTurnOff) = 0; eGEM.lb(posToTurnOff) = 0;
%% 
% Assign values of parameters for the iMAT algorithm.

% Set up hyperparameters for the linear programming section
hyperparams.eps = []; hyperparams.isgenes = true;
hyperparams.kap = []; hyperparams.eps2 = [];
hyperparams.rho = []; hyperparams.pfba = false;
hyperparams.kap2 = [];
%% 
% Set the activity coefficient for demand reactions to be a fixed value, based 
% on the medium perturbation experiments I ran with this current iteration of 
% the eGEMM.

% Get list of demand reactions from eGEMM
%pos = [3773:3777]; % Bulk reactions only
pos = [3773:3777, 3780:3799];
dm_reactions = string(eGEM.rxns(pos))';
% Activity coefficient for demand reactions
epsilon_methylation = 0.069;
epsilon_acetylation = 0.069;
for i = 1:length(dm_reactions)
    if contains(lower(dm_reactions(i)), "ac")
        epsilon(i, 1) = epsilon_acetylation;
    else
        epsilon(i, 1) = epsilon_methylation;
    end
end
%% 3. Compute base metabolic model fluxes and growth rates for comparison
% Compute the solution for an unconstrained metabolic model. This will be the 
% control model.

% for i = 1:length(dm_reactions)
%     [~, rxnPos] = ismember(dm_reactions(i), eGEM.rxns);
%     
%     % Add activity to demand reactions
%     tmp = eGEM;
%     tmp.c(rxnPos) = epsilon(i);
%     
%     [~, soln] = CFR(tmp, hyperparams, {}, {});
%     
%     grate(1, i) = soln.x(3743);
% 
%     tmp2 = tmp;
%     tmp2.lb(3743) = soln.x(3743) * 0.99;
%     tmp2.c(3743)  = 0;                                                
%     
%     [~, soln] = CFR(tmp, hyperparams, {}, {});
%     flux(1, i)  = soln.x(rxnPos);
% end
%% 4. Preprocess Cancer Cell Line Encyclopedia data
% Load data from the Cancer Cell Line Encyclopedia
% This module uses RNASeq data from the CCLE and performs some basic preprocessing. 
% The data was normalized to TPM. There are 864 cancer cell lines by 1020 metabolic 
% genes that map onto RECON1 using Entrez identifiers.

load /nfs/turbo/umms-csriram/scampit/Data/RNASeq/CCLE/CCLE_RNASeq.mat
load /nfs/turbo/umms-csriram/scampit/Data/Proteomics/CCLE/CCLE_Proteomics.mat
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
%% 5. Compute CCLE metabolic fluxes from the eGEMM

%fileName = '~/Data/CBM/eGEM/CCLE_fluxes_all.mat';
%fileName = '/nfs/turbo/umms-csriram/scampit/Data/CBM/eGEM/CCLE_fluxes_writers.mat';
fileName = '/nfs/turbo/umms-csriram/scampit/Data/CBM/eGEM/CCLE_fluxes_all.mat';

%save(fileName, '-v7.3')
% Ensure that histone readers and writers are being used in the metabolic model
% First, let's get gene identifiers that are associated with histone writers 
% and erasers.

% Load up curated histone reaction map
histoneReactionFile = '/nfs/turbo/umms-csriram/scampit/Data/Reconstructions/eGEM/Epigenome-Scale Metabolic Reconstruction Map.xlsx';
histoneReactionMap = xlsread(histoneReactionFile, 'Genes');

% Edit Entrez 
histone_entrez = histoneReactionMap(:, 1);
histone_bigg = string(histone_entrez + .1);
histone_entrez = string(histone_entrez);
% SANITY CHECK: Intersecting Entrez IDs between metabolic model and dataset
% This code block just checks and sees if there are any intersecting genes with 
% the metabolic model and the RNASeq data.

% tmp1 = [];
% tmp2 = [];
% for j = 1:size(log2fc, 2)
%         
%     % Get differentially expressed genes for each cancer cell line.
%     DE{j}.name    = cell_names(j);
%     DE{j}.medium  = medium(j);
%     DE{j}.culture = cultures(j);
%     DE{j}.tissue  = tissues(j);
%     DE{j}.upreg   = entrez_ids(log2fc(:, j) > 0 & pvalue(:, j) <= 0.05);
%     DE{j}.downreg = entrez_ids(log2fc(:, j) < 0 & pvalue(:, j) <= 0.05);
%     
%     gene_list = [DE{j}.upreg; DE{j}.downreg];
%     [~, ia, ib] = intersect(string(gene_list), string(histone_entrez));
%     
%     tmp1 = [tmp1; ia];
%     tmp2 = [tmp2; ib];
%     
% end
% Compute metabolic flux profiles
% Now let's compute the data using the model. The data will automatically be 
% saved.

% Demand reactions
flux = cell(length(dm_reactions), 1);
grate = cell(length(dm_reactions), 1);
save(fileName, 'flux', 'grate');

medium_file = '/nfs/turbo/umms-csriram/scampit/Data/Medium/FINAL_MEDIUM_MAP.xlsx';

eGEM.genes = eGEM.geneEntrezID;
entrez_ids = string(entrez_ids);

hyperparams.hscore = false;
hyperparams.gamma = [];

for i = 1:length(dm_reactions)
    [~, rxnPos] = ismember(dm_reactions(i), eGEM.rxns);
    tmp = eGEM;
    tmp.c(rxnPos) = epsilon(i);
    
    % Cell lines
    for j = 1:size(log2fc, 2)
        
        % Get differentially expressed genes for each cancer cell line.
        DE{j}.name    = cell_names(j);
        DE{j}.medium  = medium(j);
        DE{j}.culture = cultures(j);
        DE{j}.tissue  = tissues(j);
        DE{j}.upreg   = entrez_ids(log2fc(:, j) > 0 & pvalue(:, j) <= 0.05);
        DE{j}.downreg = entrez_ids(log2fc(:, j) < 0 & pvalue(:, j) <= 0.05);
        
        % Constrain the metabolic model to be medium-specific
        mediumModel = addMediumConstraints(tmp, medium(j), medium_file);
        
        try
            % Compute flux using iMAT 
            [mdl, soln] = CFR(mediumModel, hyperparams, ...
                              DE{j}.upreg, DE{j}.downreg);
                                           
            % Get the growth rate
            grates(j, 1) = soln.x(3743);
            
            % Ensure there's growth
            tmp2 = mediumModel;
            tmp2.lb(3743) = soln.x(3743) * 0.99;
            tmp2.c(3743)  = 0;                                                
            
            % Compute metabolic flux when optimizing histone markers
            [~, soln] = CFR(tmp2, hyperparams, DE{j}.upreg, DE{j}.downreg);
    
            % Store data here
            fluxes(:, j) = soln.x(1:length(eGEM.rxns));
        catch ME
            grates(j, 1) = NaN;
            fluxes(:, j) = NaN(length(eGEM.rxns), 1);
        end
        
    end
    flux{i, 1} = fluxes;
    grate{i, 1} = grates;
    save(fileName, 'flux', 'grate', 'j', 'i', '-append');
    
end
%% 6. Writers / Erasers Ratio Parameter
% For this scheme, I need to do the following steps: 
%% 
% # Separate the RNASeq data out by whether or not the differentially expressed 
% genes are writers or erasers and ensure the algorithm is distinguishing between 
% methylation and acetylation
% # Create a scaling factor for the upper bound of the demand reaction. The 
% scaling factor $\gamma$ I settled for is the following:$\;\gamma =\frac{\frac{\textrm{counts}\left(Z_{\textrm{up},\textrm{writer}} 
% \right)}{\textrm{counts}\left(Z_{\textrm{up},\textrm{eraser}} \right)}}{\frac{\textrm{counts}\left(Z_{\textrm{down},\textrm{writer}} 
% \right)}{\textrm{counts}\left(Z_{\textrm{down},\textrm{eraser}} \right)}}$
% # For the $i^{\textrm{th}}$demand reaction, constrain its upper bound using 
% $\gamma$: |model.ub(i) = gamma*model.ub(i)|
% *Case examples*
% *Case 1: All writers are up*
% 
% $$\gamma =\frac{\frac{10}{1}}{\frac{1}{1}}=10$$
% 
% *Case 2: All erasers are up* 
% 
% $$\gamma =\frac{\frac{1}{10}}{\frac{1}{1}}=0\ldotp 1$$
% 
% *Case 3: Writers are up and erasers are down*
% 
% $$\gamma =\frac{\frac{10}{1}}{\frac{1}{10}}=100$$
% 
% *Case 4: Erasers are up and writers are down* 
% 
% $$\gamma =\frac{\frac{1}{10}}{\frac{10}{1}}=0\ldotp 01$$
% Compute metabolic flux profiles
% Now let's compute the data using the model. The data will automatically be 
% saved.

% Demand reactions
flux = cell(length(dm_reactions), 1);
grate = cell(length(dm_reactions), 1);

%fileName = '/nfs/turbo/umms-csriram/scampit/Data/CBM/eGEM/CCLE_fluxes_writers_gamma.mat';
fileName = '/nfs/turbo/umms-csriram/scampit/Data/CBM/eGEM/CCLE_fluxes_all_gamma.mat';
save(fileName, 'flux', 'grate');

eGEM.genes = eGEM.geneEntrezID;
entrez_ids = string(entrez_ids);
medium_file = '/nfs/turbo/umms-csriram/scampit/Data/Medium/FINAL_MEDIUM_MAP.xlsx';

for i = 1:length(dm_reactions)
    
    % Set the objective coefficient for the model to be 0.069.
    [~, rxnPos] = ismember(dm_reactions(i), eGEM.rxns);
    tmp = eGEM;
    tmp.c(rxnPos) = epsilon(i);
    
    for j = 1:size(log2fc, 2)
        
        % Get differentially expressed genes for each cancer cell line.
        DE{j}.name    = cell_names(j);
        DE{j}.medium  = medium(j);
        DE{j}.culture = cultures(j);
        DE{j}.tissue  = tissues(j);
        DE{j}.upreg   = entrez_ids(log2fc(:, j) > 0 & pvalue(:, j) <= 0.05);
        DE{j}.downreg = entrez_ids(log2fc(:, j) < 0 & pvalue(:, j) <= 0.05);
        
        gamma = 0;
                
        % Compute a score that reflects the proportion of up and
        % downregulated writers / erasers
        if contains(lower(dm_reactions(i)), "ac")
            writers = histone_entrez(53:61);
            erasers = histone_entrez(62:80);
        elseif contains(lower(dm_reactions(i)), "me")
            writers = histone_entrez(21:52);
            erasers = histone_entrez(1:20);
        end
        
        % Compute upregulated section / numerator
        z_up_write = sum(ismember(string(DE{j}.upreg), string(writers)), 'all') + 0.1;
        z_up_erase = sum(ismember(string(DE{j}.upreg), string(erasers)), 'all') + 0.1;
        
        % Compute downregulated section / denominator
        z_down_write = sum(ismember(string(DE{j}.downreg), string(writers))) + 0.1;
        z_down_erase = sum(ismember(string(DE{j}.downreg), string(erasers)), 'all') + 0.1;
        
        % Compute gamma using numerator / denominator
        numerator = (z_up_write / z_up_erase);
        denominator = (z_down_write / z_down_erase);
        gamma = log(numerator / denominator);
        hyperparams.gamma = gamma;
        
        % Constrain the metabolic model to be medium-specific
        mediumModel = addMediumConstraints(tmp, medium(j), medium_file);
        
        try
            % Compute flux using iMAT 
            hyperparams.hscore = false;
            [mdl, soln] = CFR(mediumModel, hyperparams, ...
                              DE{j}.upreg, DE{j}.downreg);
                                           
            % Get the growth rate
            grates(j, 1) = soln.x(3743);
            
            % Ensure there's growth
            tmp2 = mediumModel;
            tmp2.lb(3743) = soln.x(3743) * 0.99;
            tmp2.c(3743)  = 0;     
            
            hyperparams.hscore = true;
            % Compute metabolic flux when optimizing histone markers
            [~, soln] = CFR(tmp2, hyperparams, ...
                            DE{j}.upreg, DE{j}.downreg);
    
            % Store data here
            fluxes(:, j) = soln.x(1:length(eGEM.rxns));
        catch ME
            grates(j, 1) = NaN;
            fluxes(:, j) = NaN(length(eGEM.rxns), 1);
        end
    end
    flux{i, 1} = fluxes;
    grate{i, 1} = grates;
    save(fileName, 'flux', 'grate', 'j', 'i', '-append'); 
end