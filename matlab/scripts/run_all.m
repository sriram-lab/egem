%% Code to run modules
initCobraToolbox;
changeCobraSolver('gurobi');

load ./../models/eGEM.mat % minimal eGEM model 
%load ./../models/acetyl2.mat % new acetylation model

%load ./../models/recon1
%model = metabolicmodel;

%load supplementary_software_code acetylation_model
%model = acetylation_model; %Shen et al., 2019

%% Correlation values between histone markers and metabolic flux
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
%for n = 1:length(epsilon2)
%    [correl, pval] = histone_corr(model, compartment,...
%        mode, epsilon, epsilon2(n), rho, kappa, minfluxflag);
%end

% Run with the optimized objective coefficients
epsilon2_excess = [1E-6, 1E-6, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-6, 1E-6,...
    1E-5, 1E-6, 1, 1E-6, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-5];

epsilon2_depletion = [1E-6, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-5,...
    1E-6, 1E-4, 1, 1E-5, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-6];

epsilon2 = [epsilon2_excess; epsilon2_depletion];
epsilon2 = epsilon2';
[correl, pval] = histone_corr(model, compartment,...
        mode, epsilon, epsilon2(n), rho, kappa, minfluxflag);


%% Heatmap of metabolic reactions vs excess/depletion of medium coponents

% Use params for testing 
%compartment = 'n';
%epsilon2 = 1E-3;
%scaling = [];
%[excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
%    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, 'n',...
%    epsilon2, []);

% Sample some epsilon values to determine largest dynamic range
%epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];

% Use epsilon values that gave the largest dynamic range in metabolic
% fluxes 
epsilon2_excess = [1E-6, 1E-6, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-6, 1E-6,...
    1E-5, 1E-6, 1, 1E-6, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-5];

epsilon2_depletion = [1E-6, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-6, 1E-6, 1E-5,...
    1E-6, 1E-4, 1, 1E-5, 1E-6, 1, 1E-6, 1E-6, 1E-6, 1E-6, 1E-6];

epsilon2 = [epsilon2_excess; epsilon2_depletion];
epsilon2 = epsilon2';

% Run metabolic_sensitivity
    % for the switch case:
        % dyn - maximize using dynamic range of metabolic fluxes
        % grate - fix biomass to max value and optimize for histone
        % markers

% For single epsilon determination
%compartment = ['n', 'c', 'm'];
%for n = 1:length(epsilon2)
    %for m = 1:length(compartment)
%    [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
%    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, 'n',...
%    epsilon2(:, n), [], 'dyn');
    %end
%end

% For optimizing with multiple epsilon values:
[excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, 'n',...
    epsilon2, [], 'dyn');

%% Density plot
A = densityplot('eGEMn');

[x,y,z] = meshgrid(1:50, 1:20, 1:6);
for i=1:6
    surf(x(:,1,1), y(1,:,1), A(:,:,i));
    hold on;
    colorbar
end