%% Run analyses for epigenome-scale metabolic models
initCobraToolbox;
changeCobraSolver('gurobi');

%% load different models

% Minial eGEM -> does not contain other one carbon reactions
load ./../models/min.mat % minimal eGEM model
%model = eGEM_min;
% Human metabolic reconstruction 1
%load ./../models/recon1
%model = metabolicmodel;

% Current eGEM
%load ./../models/eGEM.mat

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

% SRA
% Single test case
epsilon2 = 1E-3;
exp = 'single-reaction-analysis';
[sra, ~, ~, ~, ~] = metabolic_sensitivity(min, 'n', 1E-6, 'zscore', 'sra',...
    'RPMI', []);

% Get the optimal s
medium_of_interest = {'RPMI', 'DMEM', 'L15', 'McCoy 5A', 'Iscove'};
epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    for n = 1:length(epsilon2)
        % Run all
        str =  strcat("[sra", string(n), '_', medium_of_interest(med), ", ~, ~, ~, ~] = metabolic_sensitivity(min, 'n', epsilon2(n), 'zscore', 'sra',medium_of_interest(med), []);");
        eval(str)
        % Plot all
        str = strcat("plot_heatmap(sra", string(n), '_', medium_of_interest(med), ",'sra', epsilon2(n), medium_of_interest(med))");
        eval(str)
    end
end

% Calculate epsilon2 values to use for fba/fva
% DMEM
epsilon2_excess = dynamic_range(sra1_DMEM, sra2_DMEM, sra3_DMEM, sra4_DMEM,...
    sra5_DMEM, sra6_DMEM, sra7_DMEM, "excess");
epsilon2_depletion = dynamic_range(sra1_DMEM, sra2_DMEM, sra3_DMEM, sra4_DMEM,...
    sra5_DMEM, sra6_DMEM, sra7_DMEM, "depletion");
epsilon2_dmem = [epsilon2_excess, epsilon2_depletion];

% RPMI
epsilon2_excess = dynamic_range(sra1_RPMI, sra2_RPMI, sra3_RPMI, sra4_RPMI,...
    sra5_RPMI, sra6_RPMI, sra7_RPMI, "excess");
epsilon2_depletion = dynamic_range(sra1_RPMI, sra2_RPMI, sra3_RPMI, sra4_RPMI,...
    sra5_RPMI, sra6_RPMI, sra7_RPMI, "depletion");
epsilon2_rpmi = [epsilon2_excess, epsilon2_depletion];

% L15
epsilon2_excess = dynamic_range(sra1_L15, sra2_L15, sra3_L15, sra4_L15,...
    sra5_L15, sra6_L15, sra7_L15, "excess");
epsilon2_depletion = dynamic_range(sra1_L15, sra2_L15, sra3_L15, sra4_L15,...
    sra5_L15, sra6_L15, sra7_L15, "depletion");
epsilon2_l15 = [epsilon2_excess, epsilon2_depletion];

% Run FBA with different epsilon values obtained by maximizing the dynamic range:
[~, fba, ~, excess_soln, depletion_soln] = metabolic_sensitivity(min, 'n',...
        epsilon2_rpmi, 'zscore', 'fba', 'RPMI');
[~, fba, ~, excess_soln, depletion_soln] = metabolic_sensitivity(min, 'n',...
        epsilon2_dmem, 'zscore', 'fba', 'DMEM');
[~, fba, ~, excess_soln, depletion_soln] = metabolic_sensitivity(min, 'n',...
        epsilon2_l15, 'zscore', 'fba', 'L15');

medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all
    str =  strcat("[~, fba_", lower(medium_of_interest(med)),", ~, ~, ~] =metabolic_sensitivity(min, 'n', epsilon2_", lower(medium_of_interest(med)), ", 'zscore', 'fba', medium_of_interest(med), []);");
    eval(str)
    % Plot all
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)), ",'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end
    
    
% Run FBA with all the reactions maximized to 1:
proxy = ones(20,2);
for med = 1:length(medium_of_interest)
    disp(medium_of_interest(med))
    % Run all
    str =  strcat("[~, fba_", lower(medium_of_interest(med)),", ~, ~, ~] = metabolic_sensitivity(min, 'n', proxy, 'zscore', 'fba', medium_of_interest(med), []);");
    eval(str)
    % Plot all
    str = strcat("plot_heatmap(fba_", lower(medium_of_interest(med)), ",'fba', epsilon2, medium_of_interest(med))");
    eval(str)
end

% For optimizing multiple reactions simultaneously using FBA with different
% media:
[~, media] = xlsfinfo('./../../data/uptake.xlsx');
for fil=1:length(sheets)
    [~, fba, ~, excess_soln, depletion_soln] = metabolic_sensitivity(min, 'n',...
        epsilon2, 'zscore', 'fba', media(fil));
end

% For optimizing multiple reactions simultaneously using FVA with different
% media:
%[~, media] = xlsfinfo('./../../data/uptake.xlsx');
medium_of_interest = {'RPMI', 'DMEM', 'L15'};
for med=1:length(medium_of_interest)
    % Run all
    str =  strcat("[~, ~,  fva_", lower(medium_of_interest(med)), " ~, ~] = metabolic_sensitivity(min, 'n', epsilon2_", lower(medium_of_interest(med)), ", 'zscore', 'fva', medium_of_interest(med), 99);");
    eval(str)
    % Plot all
    str = strcat("plot_heatmap(fva_", lower(medium_of_interest(med)), ",'fva', epsilon2, medium_of_interest(med))");
    eval(str)
end

%% Correlation values between histone markers and metabolic flux
% INPUTS:
    % h3marks: list of H3 marks from CCLE data (column values)
    % h3names: list of CCLE cell lines (row values)
    % h3vals: matrix containing values corresponding to h3marks and h3names

% Use params for testing
% epsilon2 = 1E-3

% Initialize params for iMAT algorithm
compartment = 'n';
mode = 1;
epsilon = 1E-3;
rho = 1;
kappa = 1E-3;
minfluxflag = 0;

% Run with multiple objective coefficients to obtain dynamic range of
% values
%epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];

% Run with the optimized objective coefficients
epsilon2_excess = [1E-6, 1E-6, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-6, 1E-6,...
    1E-5, 1E-6, 1, 1E-6, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-5];

epsilon2_depletion = [1E-6, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-5,...
    1E-6, 1E-4, 1, 1E-5, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-6];

epsilon2 = [epsilon2_excess; epsilon2_depletion];
epsilon2 = epsilon2';
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