%% Run analyses for epigenome-scale metabolic models
initCobraToolbox;
changeCobraSolver('gurobi');
%% Load genome-scale metabolic models. Must run from `/matlab/scripts` directory

% The bulk eGEM model
load ./../models/eGEM.mat 
model = eGEM;

% New acetylation model containing ACSS2 and PDH in nucleus
%load ./../models/acetyl2.mat

% Minial eGEM that does not contain one carbon reactions and demethylation
% rxns
%load ./../models/min.mat % minimal eGEM model
%model = eGEM_min;

% Human metabolic reconstruction 1 (RECON1; Duarte et al., 2007)
%load ./../models/recon1
%model = metabolicmodel;

% Acetylation model from Shen et al., 2019
%load ./../shen-et-al/supplementary_software_code acetylation_model
%model = acetylation_model; 

%% Metabolic sensitivity analysis for excess and depleted medium components
    % INPUT:
        % switch case arguments:
            % single-reaction-analysis - optimizes the flux of a single
            % reaction
            % dyn - optimizes the flux of several reactions using the results
            % from `single-reaction-analysis`
            % grate - fix biomass to max value and optimize for histone
            % markers
    % OUTPUT: 
        % Heatmaps of the metabolic fluxes, shadow prices, and reduced
        % costs corresponding to each reaction / medium component pair.

% All reactions
load('./../vars/metabolites.mat')
        
% Optimization 1A: Run Single reaction activity (SRA)
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    for n = 1:length(epsilon2)
        % Run all
        str =  strcat("[sra", string(n), '_', medium_of_interest(med),...
            "] = metabolic_sensitivity(model, metabolites, 'n',", ...
            "epsilon2(n), 'zscore', 'sra', medium_of_interest(med), []);");
        eval(str)
        % Plot all
        str = strcat("plot_heatmap(sra", string(n), '_',...
            medium_of_interest(med), ", metabolites'sra', epsilon2(n), medium_of_interest(med))");
        eval(str)
    end
end

% Calculate epsilon2 values to use for fba by dynamic range
epsilon2_dmem = dynamic_range(sra1_DMEM, sra2_DMEM, sra3_DMEM, sra4_DMEM,...
    sra5_DMEM, sra6_DMEM, sra7_DMEM, "dynamic");
epsilon2_rpmi = dynamic_range(sra1_RPMI, sra2_RPMI, sra3_RPMI, sra4_RPMI,...
    sra5_RPMI, sra6_RPMI, sra7_RPMI, "dynamic");
epsilon2_l15 = dynamic_range(sra1_L15, sra2_L15, sra3_L15, sra4_L15,...
    sra5_L15, sra6_L15, sra7_L15, "dynamic");

% Optimization 2A: Run Flux balance analysis (FBA) w/ and w/o competition
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all without competition for all reactions
    str =  strcat("[fba_", lower(medium_of_interest(med)),"_noComp]", ...
        "= metabolic_sensitivity(model, metabolites, 'n', epsilon2_", ...
        lower(medium_of_interest(med)), ", 'zscore', 'no_competition',", ...
        "medium_of_interest(med), []);");
    eval(str)

    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)),...
        "_noComp, metabolites, 'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all w/ competition for all reactions
    str =  strcat("[fba_", lower(medium_of_interest(med)),"_comp, ", ...
        "] = metabolic_sensitivity(model, metabolites, 'n', epsilon2_", ...
        lower(medium_of_interest(med)), ", 'zscore', 'competition',", ...
        "medium_of_interest(med), []);"); 
    eval(str)
    
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)), ...
        "_comp, metabolites, 'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

% Optimization 3A: Run Flux variability analysis (FVA) for all reactions
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med=1:length(medium_of_interest)
    % Run FVA for all reactions
    str =  strcat("[fva_", lower(medium_of_interest(med)),...
        "] = metabolic_sensitivity(model, metabolites, 'n', epsilon2_",...
        lower(medium_of_interest(med)), ", 'zscore', 'fva', medium_of_interest(med), 90);");
    eval(str)
    
    % Plot all reactions
    str = strcat("plot_heatmap(fva_", lower(medium_of_interest(med)),...
        ", metabolites, 'fva', epsilon2, medium_of_interest(med))");
    eval(str)
end

% Save results for all reactions
[fba_nocomp_rpmi_excess, fba_nocomp_rpmi_depletion] = metabolite_dict(fba_rpmi_noComp, 'RPMI', 'T2| All Rxns FBA', 'no_competition');
[fba_nocomp_dmem_excess, fba_nocomp_dmem_depletion] = metabolite_dict(fba_dmem_noComp, 'DMEM', 'T2| All Rxns FBA', 'no_competition');
[fba_nocomp_l15_excess, fba_nocomp_l15_depletion] = metabolite_dict(fba_l15_noComp, 'L15', 'T2| All Rxns FBA', 'no_competition');
[fba_comp_rpmi_excess, fba_comp_rpmi_depletion] = metabolite_dict(fba_rpmi_comp, 'RPMI', 'T2| All Rxns FBA', 'competition');
[fba_comp_dmem_excess, fba_comp_dmem_depletion] = metabolite_dict(fba_dmem_comp, 'DMEM', 'T2| All Rxns FBA', 'competition');
[fba_comp_l15_excess, fba_comp_l15_depletion] = metabolite_dict(fba_l15_comp, 'L15', 'T2| All Rxns FBA', 'competition');
[fba_fva_rpmi_excess, fba_fva_rpmi_depletion] = metabolite_dict(fva_rpmi, 'RPMI', 'T3| All Rxns FVA', 'fva');
[fba_fva_dmem_excess, fba_fva_dmem_depletion] = metabolite_dict(fva_dmem, 'DMEM', 'T3| All Rxns FVA', 'fva');
[fba_fva_l15_excess, fba_fva_l15_depletion] = metabolite_dict(fva_l15, 'L15', 'T3| All Rxns FVA', 'fva');

% Histone reactions only
histone_rxns_only = metabolites(2:5, :);

% Optimization 1B: Run Single reaction activity (SRA) for histone reactions
% only
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    for n = 1:length(epsilon2)
        % Run all
        str =  strcat("[sra", string(n), '_', medium_of_interest(med),...
            "] = metabolic_sensitivity(model, histone_rxns_only, 'n',", ...
            "epsilon2(n), 'zscore', 'sra', medium_of_interest(med), []);");
        eval(str)
        % Plot all
        str = strcat("plot_heatmap(sra", string(n), '_',...
            medium_of_interest(med), ", histone_rxns_only, 'sra', epsilon2(n), medium_of_interest(med))");
        eval(str)
    end
end

% Calculate epsilon2 values to use for fba by dynamic range
epsilon2_dmem = dynamic_range(sra1_DMEM, sra2_DMEM, sra3_DMEM, sra4_DMEM,...
    sra5_DMEM, sra6_DMEM, sra7_DMEM, "dynamic");
epsilon2_rpmi = dynamic_range(sra1_RPMI, sra2_RPMI, sra3_RPMI, sra4_RPMI,...
    sra5_RPMI, sra6_RPMI, sra7_RPMI, "dynamic");
epsilon2_l15 = dynamic_range(sra1_L15, sra2_L15, sra3_L15, sra4_L15,...
    sra5_L15, sra6_L15, sra7_L15, "dynamic");

% Optimization 2B: Run Flux balance analysis (FBA) w/ and w/o competition
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    
    % Run all without competition for histone reactions
    str =  strcat("[fba_", lower(medium_of_interest(med)),"_noComp]", ...
        "= metabolic_sensitivity(model, histone_rxns_only, 'n', epsilon2_", ...
        lower(medium_of_interest(med)), ", 'zscore', 'no_competition',", ...
        "medium_of_interest(med), []);");
    eval(str)
    
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)),...
        "_noComp, histone_rxns_only, 'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all w/ competition for histone reactions
    str =  strcat("[fba_", lower(medium_of_interest(med)),"_comp", ...
        "] = metabolic_sensitivity(model, histone_rxns_only, 'n', epsilon2_", ...
        lower(medium_of_interest(med)), ", 'zscore', 'competition',", ...
        "medium_of_interest(med), []);"); 
    eval(str)
    
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)), ...
        "_comp, histone_rxns_only, 'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

% Optimization 3B: Run Flux variability analysis (FVA) for histone
% reactions
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med=1:length(medium_of_interest)
    str =  strcat("[fva_", lower(medium_of_interest(med)),...
        "] = metabolic_sensitivity(model, histone_rxns_only, 'n', epsilon2_",...
        lower(medium_of_interest(med)), ", 'zscore', 'fva', medium_of_interest(med), 100);");
    eval(str)
    
    % Plot
    str = strcat("plot_heatmap(fva_", lower(medium_of_interest(med)),...
        ", histone_rxns_only, 'fva', epsilon2, medium_of_interest(med))");
    eval(str)
end

% Save results for histone reactions only
[fba_nocomp_rpmi_excess, fba_nocomp_rpmi_depletion] = metabolite_dict(fba_rpmi_noComp, histone_rxns_only, 'RPMI', 'T2| All Rxns FBA', 'no_competition');
[fba_nocomp_dmem_excess, fba_nocomp_dmem_depletion] = metabolite_dict(fba_dmem_noComp, histone_rxns_only, 'DMEM', 'T2| All Rxns FBA', 'no_competition');
[fba_nocomp_l15_excess, fba_nocomp_l15_depletion] = metabolite_dict(fba_l15_noComp, histone_rxns_only, 'L15', 'T2| All Rxns FBA', 'no_competition');
[fba_comp_rpmi_excess, fba_comp_rpmi_depletion] = metabolite_dict(fba_rpmi_comp, histone_rxns_only, 'RPMI', 'T2| All Rxns FBA', 'competition');
[fba_comp_dmem_excess, fba_comp_dmem_depletion] = metabolite_dict(fba_dmem_comp, histone_rxns_only, 'DMEM', 'T2| All Rxns FBA', 'competition');
[fba_comp_l15_excess, fba_comp_l15_depletion] = metabolite_dict(fba_l15_comp, histone_rxns_only, 'L15', 'T2| All Rxns FBA', 'competition');
[fba_fva_rpmi_excess, fba_fva_rpmi_depletion] = metabolite_dict(fva_rpmi, histone_rxns_only, 'RPMI', 'T3| All Rxns FVA', 'fva');
[fba_fva_dmem_excess, fba_fva_dmem_depletion] = metabolite_dict(fva_dmem, histone_rxns_only, 'DMEM', 'T3| All Rxns FVA', 'fva');
[fba_fva_l15_excess, fba_fva_l15_depletion] = metabolite_dict(fva_l15, histone_rxns_only, 'L15', 'T3| All Rxns FVA', 'fva');

%% Correlation values between histone markers and metabolic flux
% INPUTS:
    % h3marks: list of H3 marks from CCLE data (column values)
    % h3names: list of CCLE cell lines (row values)
    % h3vals: matrix containing values corresponding to h3marks and h3names

% Initialize params for iMAT algorithm
compartment = 'n';
mode = 1;
epsilon = 1E-3;
rho = 1;
kappa = 1E-3;
minfluxflag = 0;
    
[correl, pval] = histone_corr(model, compartment,...
        mode, epsilon, epsilon2(n), rho, kappa, minfluxflag);
    
%% Density plot
% A = densityplot('eGEMn');
% 
% [x,y,z] = meshgrid(1:50, 1:20, 1:6);
% for i=1:6
%     surf(x(:,1,1), y(1,:,1), A(:,:,i));
%     hold on;
%     colorbar
% end