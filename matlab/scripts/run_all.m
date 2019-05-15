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

h3marks = textread('C:\Users\scampit\Desktop\MeGEM\matlab\data\h3marks.txt',...
    '%s', 'whitespace', '\n');
h3names = textread('C:\Users\scampit\Desktop\MeGEM\matlab\data\ccle_names.txt',...
    '%s', 'whitespace', '\n');
h3vals = dlmread('C:\Users\scampit\Desktop\MeGEM\matlab\data\h3_relval.txt',...
    ',');

histone_corr(model, 'amet', [], 'n', 1, 1E-2, 1, 1E-3, 0);
rxnpos  = [find(ismember(model.rxns, 'EX_KAC'));];

%% Heatmap of metabolic reactions vs excess/depletion of medium coponents

% Use params for testing 
%compartment = 'n';
%epsilon2 = 1E-3;
%scaling = [];
%[excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
%    excess_shadow, depletion_shadow] = make_heatmap(model, 'n',...
%    epsilon2, []);

epsilon2 = [1E-4, 1E-3, 1E-2, 0.1, 1];
%compartment = ['n', 'c', 'm'];
for n = 1:length(epsilon2)
    %for m = 1:length(compartment)
    [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = make_heatmap(model, 'c',...
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