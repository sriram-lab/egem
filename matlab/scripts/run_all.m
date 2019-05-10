%% Code to run modules
initCobraToolbox;
changeCobraSolver('gurobi');

% load different models
load ./../models/model.mat % minimal eGEM model - need to 

%load recon1
%model = metabolicmodel;

%load supplementary_software_code acetylation_model
%model = acetylation_model; %Shen et al., 2019

%% Correlation values between histone markers and metabolic flux
%histone_corr(model, 'amet', [], 'n', 1, 1E-2, 1, 1E-3, 0);
%rxnpos  = [find(ismember(model.rxns, 'EX_KAC'));];

%% Heatmap of metabolic reactions vs excess/depletion of medium coponents

% Use params for testing 
%compartment = 'n';
%epsilon2 = 1E-3;
%scaling = [];

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
end