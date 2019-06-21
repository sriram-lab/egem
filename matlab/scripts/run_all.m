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

% Use params for testing 
compartment = 'n';
epsilon2 = 1E-3;
scaling = [];

var = {'./../var/metabolites.mat', './../var/cellmedia.mat', './../var/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk})
end

% Input filename for saving --> !Add this as argument to function later!
filename = './../tables/eGEM.xlsx';
colname = metabolites(:,3)';
rowname = mediareactions1(:,2);

% Metabolic sensitivity analysis for optimized epsilon values obtained from
% single optimization steps:

% Single test case for epsilon
epsilon2 = 1E-3
exp = 'single-reaction-analysis'

epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
%compartment = ['n', 'c', 'm'];
%for n = 1:length(epsilon2)
    %for m = 1:length(compartment)
[sra, ~, ~] = metabolic_sensitivity(min, 'n',...
1E-6, [], 'sra');
    %end
%end

% Use epsilon values that gave the largest dynamic range in metabolic
% fluxes. !This was manually done, but should also be codified at some point! 
epsilon2_excess = [1E-6, 1E-6, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-6, 1E-6,...
    1E-5, 1E-6, 1, 1E-6, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-5];

epsilon2_depletion = [1E-6, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-5,...
    1E-6, 1E-4, 1, 1E-5, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-6];

epsilon2 = [epsilon2_excess; epsilon2_depletion];
epsilon2 = epsilon2';

% For optimizing multiple reactions simultaneously using FBA:
[~, fba, ~, excess_soln, depletion_soln] = metabolic_sensitivity(min, 'n',...
    epsilon2, 'zscore', 'fba', 'RPMI');

% For optimizing multiple reactions simultaneously using FVA
[~, ~, fva, ~] = metabolic_sensitivity(min, 'n',...
    epsilon2, 'zscore', 'fva', 'RPMI', 99);

%% Density plot
A = densityplot('eGEMn');

[x,y,z] = meshgrid(1:50, 1:20, 1:6);
for i=1:6
    surf(x(:,1,1), y(1,:,1), A(:,:,i));
    hold on;
    colorbar
end