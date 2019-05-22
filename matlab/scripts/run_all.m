%% Code to run modules
initCobraToolbox;
changeCobraSolver('gurobi');

load ./../models/eGEM.mat % minimal eGEM model 
load ./../models/acetyl2.mat % new acetylation model

load ./../models/recon1
model = metabolicmodel;

%load supplementary_software_code acetylation_model
%model = acetylation_model; %Shen et al., 2019

%% Correlation values between histone markers and metabolic flux
% h3marks: list of H3 marks from CCLE data (column values)
% h3names: list of CCLE cell lines (row values)
% h3vals: matrix containing values corresponding to h3marks and h3names

% initialize parameters:
compartment = 'n';
mode = 1;
epsilon = 1E-3;
rho = 1;
kappa = 1E-3;
minfluxflag = 0;
epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
for n = 1:length(epsilon2)
    [correl, pval] = histone_corr(model, compartment,...
        mode, epsilon, epsilon2(n), rho, kappa, minfluxflag);
end

%% Heatmap of metabolic reactions vs excess/depletion of medium coponents

% Use params for testing 
compartment = 'n';
epsilon2 = 1E-3;
scaling = [];
[excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, 'n',...
    epsilon2, []);

epsilon2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1, 1];
%compartment = ['n', 'c', 'm'];
for n = 1:length(epsilon2)
    %for m = 1:length(compartment)
    [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, 'n',...
    epsilon2(n), []);
    %end
end

%% Density plot
A = densityplot('eGEMn');

[x,y,z] = meshgrid(1:50, 1:20, 1:6);
for i=1:6
    surf(x(:,1,1), y(1,:,1), A(:,:,i));
    hold on;
    colorbar
end