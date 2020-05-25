%% eGEM hyperparameter optimization
%% Summary
% The goal of this module is to identify the best histone activity coefficient 
% |zeta| that maximize the variance between histone metabolic fluxes in the epigenome-scale 
% metabolic model.
% Parameter distribution
% First, let's create a random distribution of values using a log distribution 
% ranging from 0 to 1. This will be used for the following modules.

clear all;
zeta_dist = lognrnd(1E-3, 1, [1, 100000]);
zeta_dist = 10.^-zeta_dist;
%histogram(zeta_dist, 100);
%xlabel('Range of values'); ylabel('Value bins (log-scale)')
%set(gca, 'YScale', 'log'); title('Parameter distribution')
% Default model settings
% I will turn off some reactions that have been known to cause issues with the 
% flux distribution, and will make the lower bound of the growth to be 0, while 
% setting the biomass objective coefficient to 1. . 

initCobraToolbox; changeCobraSolver('gurobi', 'all');
load 05042020eGEM.mat; model = eGEM;

% Turn off reactions that negatively affect model performance.
posToTurnOff = find(ismember(model.rxns, {'ALR', 'MGSA', 'MGSA2'}));
model.ub(posToTurnOff) = 0; model.lb(posToTurnOff) = 0;
model.lb(ismember(model.rxns, 'biomass_objective')) = 0; % Change the lower bound of the biomass objective function to 0. 

% Reset the primary objective function
model.c = zeros(size(model.c));
model.c(ismember(model.rxns, 'biomass_objective')) = 1;

% Replace the primary objective function with ATP production
%model = addExchangeRxn(model, {'atp_c'}, 0, 100);
%model.c(find(model.c)) = 0; model.c(ismember(model.rxns, 'biomass_objective')) = 1;
%% Make COBRA model into Gurobi form
% Next, I'll create some variables in the COBRA structure that are needed to 
% use the gurobi optimzer.

mediumMap = 'FINAL_MEDIUM_MAP.xlsx';
reactions_to_optimize = ["DM_KAC", "DM_KMe1", "DM_KMe2", "DM_KMe3"];
model.lb(ismember(model.rxns, reactions_to_optimize)) = 0;
model.rev(ismember(model.rxns, reactions_to_optimize)) = 0;

% Set up Gurobi model
model.A   = model.S;
model.obj = model.c;
model.rhs = model.b;

if (exist('model.csense','var')) && (~isempty(model.csense))
    model.sense = model.csense;
    model.sense(ismember(model.sense,'E')) = '=';
    model.sense(ismember(model.sense,'L')) = '<';
    model.sense(ismember(model.sense,'G')) = '>';
else
    model.sense = repmat('=', [size(model.S, 1), 1]);
end
model.vtype      = repmat('C', size(model.S, 2), 1);
model.modelsense = 'max';

params.outputflag          = 0;
params.Threads             = 4;
params.Seed                = 314;
params.NumericFocus        = 3;
%% Nonspecific histone flux optimization
% This block of code compute the histone fluxes if they are optimized using 
% an arbitrary activity coefficient of 1E-3.

model = addMediumConstraints(model, 'RPMI', mediumMap);
for j = 1:length(reactions_to_optimize)
    tmp_mdl = model;
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize(j))) = 1E-3;
    soln{j} = optimizeCbModel(tmp_mdl);
    grate(j) = soln{j}.x(find(model.c));
    flux(j) = soln{j}.x(ismember(tmp_mdl.rxns, reactions_to_optimize(j)));  
end
%% Load gene expression data for the CCLE 
% This is microarray data from Shen et al., 2019. I'll make it to the global 
% chromatin profiles, and compute the |zeta| values for those intersecting cancer 
% cell lines.

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
%% Histone activity coefficients across different medium
% Single histone reaction optimization 
% This module will get the parameters that maximize the variance of histone 
% metabolic fluxes under specific nutrient conditions. We'll optimize individual 
% histone marker reaction fluxes to get these parameters. 
% Determine best |zeta| for each medium condition
% First we'll get the |zeta| value that maximizes the histone flux variance 
% for individual histone reactions.

% Single reaction optimization
[~, sheetNames] = xlsfinfo(mediumMap);                                       % Read in each medium condition
final_zeta = zeros([length(sheetNames), length(reactions_to_optimize) ]);
for i = 1:length(sheetNames)                                                 % Create individual medium models
    mediumModel = addMediumConstraints(model, sheetNames(i), mediumMap);
    
    for j =  1:length(reactions_to_optimize)                                 % Run for all 4 reactions
        
        % Compute some fluxes with varying zeta
        for k = 1:1000                                                        % Sample 1000 different psi values
            tmp_mdl = mediumModel;
            zeta(k) = datasample(zeta_dist, 1);                          
            tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize(j))) = zeta(k);
            tmp_mdl.obj(ismember(tmp_mdl.rxns, reactions_to_optimize(j))) = zeta(k);
            soln = gurobi(tmp_mdl, params);                              
            histone_flux(k, 1) = soln.x(ismember(tmp_mdl.rxns, reactions_to_optimize(j)));  % Growth rates and histone fluxes
        end
        
        % Filter values that get you fluxes you want
        no_growth = (histone_flux == 0);                               % No growth
        histone_flux(no_growth) = []; zeta(no_growth) = [];
        
        no_flux = (histone_flux == 0);                                 % No flux in that specific reaction
        histone_flux(no_flux) = []; zeta(no_flux) = [];
                
        flux_var = var(histone_flux(:, 2:end), [], 2);                       % Without grate, get idx with max variance in histone fluxes      
        [~, idx] = max(flux_var, [], 1);
        
        final_zeta(i, j) = zeta(idx);                                        % Capture the zeta value that has the max var(histone_flux)
        save('zeta.mat', 'final_zeta');
    end
end
% Save the |zeta| parameters from single reaction optimization

% Save Parameters
tbl = array2table(final_zeta, ...
    'RowNames', sheetNames, ...
    'VariableNames', reactions_to_optimize);
zeta = tbl;
save('zeta_sra.mat', 'zeta');
% Use the |zeta| values to compute a flux distribution

hyperparams.eps = [];
hyperparams.kap  = [];
hyperparams.rho  = [];
hyperparams.mode = true;
hyperparams.eps2 = [];
hyperparams.pfba = true;

for i = 1:length(ia)
    constrainedModel = model;
    mediumModel = addMediumConstraints(constrainedModel, string(gcp_medium(i)));     % Add on medium constraints
    diffExpGenes = findDiffExpGenes(mediumModel, string(gcp_cellnames(i)));          % Find differentially expressed genes for this cancer cell line
    ON_fieldname = string(strcat('ON_', string(gcp_cellnames(i))));                  % Make fieldnames for structure
    OFF_fieldname = string(strcat('OFF_', string(gcp_cellnames(i))));
    
    tmp_mdl = mediumModel;
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(string(gcp_medium(i)), :);
    [~, soln] = CFR(tmp_mdl, hyperparams, ...
                    diffExpGenes.(ON_fieldname), diffExpGenes.(OFF_fieldname));
    histone_flux(k, :) = soln.x([3743; find(ismember(tmp_mdl.rxns, reactions_to_optimize))]);          % Growth rates and histone fluxes
end
% Save the histone fluxes from single reaction optimization

tbl = array2table(histone_flux, ...
    'RowNames', gcp_cellnames, ...
    'VariableNames', ["Biomass", reactions_to_optimize]);
histone_flux = tbl;
save('ccle_flux_sra.mat', 'histone_flux');
% Combination histone reaction optimization
% This code now samples histone reaction coefficients |zeta| assuming that all 
% variables are being optimized simultaneously. The |zeta| values that have the 
% max variance among histone reactions is chosen for each medium condition.
% Determine best |zeta| for each medium condition
% First we'll get the |zeta| value that maximizes the histone flux variance 
% for individual histone reactions.

clear zeta
[~, sheetNames] = xlsfinfo(mediumMap);                                       % Read in each medium condition
final_zeta = zeros([length(sheetNames), length(reactions_to_optimize)]);
for i = 1:length(sheetNames)                                                 % Create individual medium models
    mediumModel = addMediumConstraints(model, sheetNames(i), mediumMap);
            
    % Compute some fluxes with varying zeta
    for k = 1:1000                                                        % Sample 1000 different psi values
        tmp_mdl = mediumModel;
        zeta(k, :) = datasample(zeta_dist, 4);                        % Set Zeta to be a value 
        tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(k, :);
        tmp_mdl.obj(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(k, :);
        soln = gurobi(tmp_mdl, params);                              % FBA
        histone_flux(k, :) = soln.x(ismember(tmp_mdl.rxns, reactions_to_optimize));  % Growth rates and histone fluxes
    end
    
    % Filter values that get you fluxes you want
    no_growth = (histone_flux(:, 1) == 0);                               % No growth
    histone_flux(no_growth, :) = [];
    zeta(no_growth, :) = [];
            
    flux_var = var(histone_flux(:, 2:end), [], 2);                       % Without grate, get idx with max variance in histone fluxes      
    [~, idx] = max(flux_var, [], 1);
    final_zeta(i, :) = zeta(idx, :);                                        % Capture the zeta value that has the max var(histone_flux)
    save('zeta_cra.mat', 'final_zeta'); 
end
% Save the histone fluxes from competition

tbl = array2table(final_zeta, ...
    'RowNames', sheetNames, ...
    'VariableNames', reactions_to_optimize);
zeta = tbl;
save('zeta_cra.mat', 'zeta');
% Use the |zeta| values to compute a flux distribution

hyperparams.eps = [];
hyperparams.kap  = [];
hyperparams.rho  = [];
hyperparams.mode = true;
hyperparams.eps2 = [];
hyperparams.pfba = true;

for i = 1:length(ia)
    constrainedModel = model;
    mediumModel = addMediumConstraints(constrainedModel, string(gcp_medium(i)));     % Add on medium constraints
    diffExpGenes = findDiffExpGenes(mediumModel, string(gcp_cellnames(i)));          % Find differentially expressed genes for this cancer cell line
    ON_fieldname = string(strcat('ON_', string(gcp_cellnames(i))));                  % Make fieldnames for structure
    OFF_fieldname = string(strcat('OFF_', string(gcp_cellnames(i))));
    
    tmp_mdl = mediumModel;
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(string(gcp_medium(i)), :);
    [~, soln] = CFR(tmp_mdl, hyperparams, ...
                    diffExpGenes.(ON_fieldname), diffExpGenes.(OFF_fieldname));
    histone_flux(k, :) = soln.x([3743; find(ismember(tmp_mdl.rxns, reactions_to_optimize))]);          % Growth rates and histone fluxes
end
% Save the histone fluxes from single reaction optimization

tbl = array2table(histone_flux, ...
    'RowNames', gcp_cellnames, ...
    'VariableNames', ["Biomass", reactions_to_optimize]);
histone_flux = tbl;
save('ccle_flux_cra.mat', 'histone_flux');
% Single Reaction Flux Variability histone reaction optimization
% Next, I will calculate the fluxes using flux variability analysis, where I 
% am optimizing the growth rate at 99%, as well as the flux through a single histone 
% reaction. I will take the |zeta| value that results in a maximum growth rate 
% and that specific histone flux. 
% Determine best |zeta| for each medium condition
% First we'll get the |zeta| value that maximizes the histone flux variance 
% for individual histone reactions.

clear zeta histone_flux
[~, sheetNames] = xlsfinfo(mediumMap);                                       % Read in each medium condition
final_zeta = zeros([length(sheetNames), length(reactions_to_optimize)]);
for i = 1:length(sheetNames)                                                 % Create individual medium models
    mediumModel = addMediumConstraints(model, sheetNames(i), mediumMap);
    for j =  1:length(reactions_to_optimize)                                 % Run for all 4 reactions
        for k = 1:1000
            tmp_mdl = mediumModel;
            zeta(k) = datasample(zeta_dist, 1);                          % Set Zeta to be a value 
            tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize(j))) = zeta(k);
            tmp_mdl.obj(ismember(tmp_mdl.rxns, reactions_to_optimize(j))) = zeta(k);
            [~, maxFlux] = fluxVariability(tmp_mdl, 99, 'max', ...
                convertStringsToChars(["biomass_objective"; reactions_to_optimize(j)]), ...
                0, false);
            histone_flux(k, 1:2) = maxFlux';  % Growth rates and histone fluxes
        end
        no_growth = (histone_flux(:, 1) == 0);                           % No growth
        histone_flux(no_growth, :) = []; zeta(no_growth) = [];
        
        no_flux = histone_flux(:, 2) == 0;                               % No flux in that specific reaction
        histone_flux(no_flux, :) = []; zeta(no_flux) = [];
        
        [~, idx] = max(histone_flux(:, 2), [], 1);
        
        final_zeta(i, j) = zeta(idx);                                    % Capture the zeta value that has the max var(histone_flux) 
        save('zeta_sfva.mat', 'final_zeta');
    end
end
% Save the |zeta| parameters

tbl = array2table(final_zeta, ...
    'RowNames', sheetNames, ...
    'VariableNames', reactions_to_optimize);
zeta = tbl;
save('zeta_sfva.mat', 'zeta');
% Use the |zeta| values to compute a flux distribution

hyperparams.eps = [];
hyperparams.kap  = [];
hyperparams.rho  = [];
hyperparams.mode = true;
hyperparams.eps2 = [];
hyperparams.pfba = true;

for i = 1:length(ia)
    constrainedModel = model;
    mediumModel = addMediumConstraints(constrainedModel, string(gcp_medium(i)));     % Add on medium constraints
    diffExpGenes = findDiffExpGenes(mediumModel, string(gcp_cellnames(i)));          % Find differentially expressed genes for this cancer cell line
    ON_fieldname = string(strcat('ON_', string(gcp_cellnames(i))));                  % Make fieldnames for structure
    OFF_fieldname = string(strcat('OFF_', string(gcp_cellnames(i))));
    
    tmp_mdl = mediumModel;
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(string(gcp_medium(i)), :);
    [~, soln] = CFR(tmp_mdl, hyperparams, ...
                    diffExpGenes.(ON_fieldname), diffExpGenes.(OFF_fieldname));
    histone_flux(k, :) = soln.x([3743; find(ismember(tmp_mdl.rxns, reactions_to_optimize))]);          % Growth rates and histone fluxes
end
% Save the histone fluxes from single reaction optimization

tbl = array2table(histone_flux, ...
    'RowNames', gcp_cellnames, ...
    'VariableNames', ["Biomass", reactions_to_optimize]);
histone_flux = tbl;
save('ccle_flux_sfva.mat', 'histone_flux');
% Competitive Flux Variability histone reaction optimization
% This module performs competition flux variability analysis, where I am optimizing 
% the growth rate at 99%, as well as the flux through all 4 histone reaction. 
% I will take the |zeta| value that results in a maximum growth rate and the maximum 
% variance across histone fluxes. 
% Determine best |zeta| for each medium condition
% First we'll get the |zeta| value that maximizes the histone flux variance 
% for individual histone reactions.

clear zeta histone_flux
[~, sheetNames] = xlsfinfo(mediumMap);                                       % Read in each medium condition
for i = 1:length(sheetNames)                                                 % Create individual medium models
    mediumModel = addMediumConstraints(model, sheetNames(i), mediumMap);
    for k = 1:1000                                                        % Sample 1000 different psi values
        tmp_mdl = mediumModel;
        zeta(k, :) = datasample(zeta_dist, 4);                        % Set Zeta to be a value 
        tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(k, :);
        tmp_mdl.obj(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(k, :);
        [~, maxFlux] = fluxVariability(tmp_mdl, 99, 'max', ...
                    convertStringsToChars(["biomass_objective", reactions_to_optimize]), ...
                    0, false);
        histone_flux(k, :) = maxFlux';  % Growth rates and histone fluxes
    end
    
    no_growth = (histone_flux(:, 1) == 0);                               % No growth
    histone_flux(no_growth, :) = [];
    zeta(no_growth, :) = [];
            
    flux_var = var(histone_flux(:, 2:end), [], 2);                       % Without grate, get idx with max variance in histone fluxes      
    [~, idx] = max(flux_var, [], 1);
    final_zeta(i, :) = zeta(idx, :);                                      % Capture the zeta value that has the max var(histone_flux)
    save('zeta_cfva.mat', 'final_zeta');
end    
% Save the |zeta| parameters

tbl = array2table(final_zeta, ...
    'RowNames', sheetNames, ...
    'VariableNames', reactions_to_optimize);
zeta = tbl;
save('zeta_cfva.mat', 'zeta');
% Use the |zeta| values to compute a flux distribution

hyperparams.eps = [];
hyperparams.kap  = [];
hyperparams.rho  = [];
hyperparams.mode = true;
hyperparams.eps2 = [];
hyperparams.pfba = true;

for i = 1:length(ia)
    constrainedModel = model;
    mediumModel = addMediumConstraints(constrainedModel, string(gcp_medium(i)));     % Add on medium constraints
    diffExpGenes = findDiffExpGenes(mediumModel, string(gcp_cellnames(i)));          % Find differentially expressed genes for this cancer cell line
    ON_fieldname = string(strcat('ON_', string(gcp_cellnames(i))));                  % Make fieldnames for structure
    OFF_fieldname = string(strcat('OFF_', string(gcp_cellnames(i))));
    
    tmp_mdl = mediumModel;
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = zeta(string(gcp_medium(i)), :);
    [~, soln] = CFR(tmp_mdl, hyperparams, ...
                    diffExpGenes.(ON_fieldname), diffExpGenes.(OFF_fieldname));
    histone_flux(k, :) = soln.x([3743; find(ismember(tmp_mdl.rxns, reactions_to_optimize))]);          % Growth rates and histone fluxes
end
% Save the histone fluxes from single reaction optimization

tbl = array2table(histone_flux, ...
    'RowNames', gcp_cellnames, ...
    'VariableNames', ["Biomass", reactions_to_optimize]);
histone_flux = tbl;
save('ccle_flux_cfva.mat', 'histone_flux');