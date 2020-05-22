%% Competition reaction activity analysis <SERVER VERSION>

%% Create a distribution of zeta values to sample
clear all;
zeta_dist = lognrnd(1E-3, 1, [1, 100000]);
zeta_dist = 10.^-zeta_dist;
%histogram(zeta_dist, 100);
%xlabel('Range of values'); ylabel('Value bins (log-scale)')
%set(gca, 'YScale', 'log'); title('Parameter distribution')

%% Load model and COBRA toolbox
initCobraToolbox; changeCobraSolver('gurobi', 'all');
load 05042020eGEM.mat; model = eGEM;

% Turn off reactions that negatively affect model performance.
posToTurnOff = find(ismember(model.rxns, {'ALR', 'MGSA', 'MGSA2'}));
model.ub(posToTurnOff) = 0; model.lb(posToTurnOff) = 0;
model.lb(ismember(model.rxns, 'biomass_objective')) = 0; % Change the lower bound of the biomass objective function to 0. 

% Reset the primary objective function
model.c = zeros(size(model.c));
model.c(ismember(model.rxns, 'biomass_objective')) = 1;


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

%% Save the zeta values
tbl = array2table(final_zeta, ...
    'RowNames', sheetNames, ...
    'VariableNames', reactions_to_optimize);
zeta = tbl;
save('zeta_cra.mat', 'zeta');

%% Compute fluxes using the zeta values
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
    tmp_mdl.c(ismember(tmp_mdl.rxns, reactions_to_optimize)) = table2array(zeta(string(gcp_medium(i)), :));
    [~, soln] = CFR(tmp_mdl, hyperparams, ...
                    diffExpGenes.(ON_fieldname), diffExpGenes.(OFF_fieldname));
    histone_flux(k, :) = soln.x([3743; find(ismember(tmp_mdl.rxns, reactions_to_optimize))]);          % Growth rates and histone fluxes
end

%% Save histone fluxes
tbl = array2table(histone_flux, ...
    'RowNames', gcp_cellnames, ...
    'VariableNames', ["Biomass", reactions_to_optimize]);
histone_flux = tbl;
save('ccle_flux_cra.mat', 'histone_flux');

