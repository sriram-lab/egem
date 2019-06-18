%% Code to run modules
initCobraToolbox;
changeCobraSolver('gurobi');

<<<<<<< HEAD
% load different models
load ./../models/model.mat % minimal eGEM model - need to 

%load recon1
=======
load ./../models/eGEM.mat % minimal eGEM model 
%load ./../models/acetyl2.mat % new acetylation model

%load ./../models/recon1
>>>>>>> ff015e3c4b46fc0ff3510c7840d25949dba527ce
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
<<<<<<< HEAD

%[excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
%    excess_shadow, depletion_shadow] = make_heatmap(model, 'n',...
%    epsilon2, []);

var = {'./../var/metabolites.mat', './../var/cellmedia.mat', './../var/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk})
end

% input filename for saving
filename = './../tables/eGEM.xlsx';
colname = metabolites(:,3)';
rowname = mediareactions1(:,2);

epsilon2 = [1E-4, 1E-3, 1E-2, 0.1, 1, 1E4];
for n=1:length(epsilon2)
    [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = make_heatmap(model, 'n',...
    epsilon2(n), []);
    % Excess flux
    xlswrite(filename, colname, string(epsilon2(n)), 'B1:U1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A2:A51');
    xlswrite(filename, excess_flux, string(epsilon2(n)), 'B2:U51');
    % Depleted flux
    xlswrite(filename, colname, string(epsilon2(n)), 'X1:AQ1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A2:A51');
    xlswrite(filename, depletion_flux, string(epsilon2(n)), 'X2:AQ51');
    % Excess reduced cost
    xlswrite(filename, colname, string(epsilon2(n)), 'B1:U1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A54:A103');
    xlswrite(filename, excess_redcost, string(epsilon2(n)), 'B54:U103');
    % Depleted reduced cost
    xlswrite(filename, colname, string(epsilon2(n)), 'X1:AQ1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A54:A103');
    xlswrite(filename, depletion_redcost, string(epsilon2(n)), 'X54:AQ103');
    % Excess shadow price
    xlswrite(filename, colname, string(epsilon2(n)), 'B1:U1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A106:A155');
    xlswrite(filename, excess_shadow, string(epsilon2(n)), 'B106:U155');
    % Depleted shadow price
    xlswrite(filename, colname, string(epsilon2(n)), 'X1:AQ1');
    xlswrite(filename, rowname, string(epsilon2(n)), 'A106:A155');
    xlswrite(filename, depletion_shadow, string(epsilon2(n)), 'X106:AQ155');
=======
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
>>>>>>> ff015e3c4b46fc0ff3510c7840d25949dba527ce
end