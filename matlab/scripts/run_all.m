%% Run analyses for epigenome-scale metabolic models
initCobraToolbox;
changeCobraSolver('gurobi');
%% load models. Must run from "scripts" directory
load ./../models/eGEM.mat % bulk eGEM model 
model = eGEM;

%load ./../models/acetyl2.mat % new acetylation model
% Minial eGEM -> does not contain other one carbon reactions
%load ./../models/min.mat % minimal eGEM model
%model = eGEM_min;
% Human metabolic reconstruction 1
%load ./../models/recon1
%model = metabolicmodel;

% New acetylation metabolic model. Contains PDC and ACSS2.
%load ./../models/acetyl2.mat % new acetylation model

% Old acetylation metabolic model
%load supplementary_software_code acetylation_model
%model = acetylation_model; %Shen et al., 2019
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

% Optimization 1: Run Single reaction activity (SRA)
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    for n = 1:length(epsilon2)
        % Run all
        str =  strcat("[sra", string(n), '_', medium_of_interest(med), ", ~, ~, ~, ~] = metabolic_sensitivity(model, 'n', epsilon2(n), 'zscore', 'sra', medium_of_interest(med), []);");
        eval(str)
        % Plot all
        str = strcat("plot_heatmap(sra", string(n), '_', medium_of_interest(med), ",'sra', epsilon2(n), medium_of_interest(med))");
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

% Optimization 2: Run Flux balance analysis (FBA)
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all
    str =  strcat("[~, fba_", lower(medium_of_interest(med)),", ~, ~, ~] =metabolic_sensitivity(model, 'n', epsilon2_", lower(medium_of_interest(med)), ", 'zscore', 'fba', medium_of_interest(med), []);");
    eval(str)
    % Plot all
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)), ",'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

% Optimization 3: Run Flux variability analysis (FVA)
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med=1:length(medium_of_interest)
    % Run all
    str =  strcat("[~, ~,  fva_", lower(medium_of_interest(med)), " ~, ~] = metabolic_sensitivity(model, 'n', epsilon2_", lower(medium_of_interest(med)), ", 'zscore', 'fva', medium_of_interest(med), 99);");
    eval(str)
    % Plot all
    str = strcat("plot_heatmap(fva_", lower(medium_of_interest(med)), ",'fva', epsilon2, medium_of_interest(med))");
    %eval(str)
end

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
A = densityplot('eGEMn');

[x,y,z] = meshgrid(1:50, 1:20, 1:6);
for i=1:6
    surf(x(:,1,1), y(1,:,1), A(:,:,i));
    hold on;
    colorbar
end