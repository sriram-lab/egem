%% histone_corr calculates the correlation value between various histone markers and the metabolic flux obtained from the iMAT algorithm
function [rho, pval] = histone_corr(model, compartment, mode, epsilon, epsilon2, rho, kappa, minfluxflag)
%% INPUTS:
    % model: Initial genome scale model
    % compartment: Subcellular compartment of interest
    % mode: for constrain flux regulation
    % epsilon: for constrain flux regulation
    % epsilon2: obj coef weight for reaction of interest
    % rho: for constrain flux regulation
    % kappa: for constrain flux regulation
    % minfluxflag: for parsimonious flux balance analysis
%% OUTPUTS:
    % correl: correlation values associated with each histone marker/rxn 
    % pval: the p-value associated with correl
    % cell_line_match: cell lines that matched between gene expression and
    % proteomics
    % heatmap that visualizes the correlation values
%% histone_corr
load ./../vars/supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz % contains CCLE cellline names for gene exp, enzymes encoding specific metabolites, and gene expression data (z-transformed)

% New variables
path = './../new_var/';
vars = {...
    [path 'h3_ccle_names.mat'], [path 'h3_marks.mat'],...
    [path 'h3_media.mat'], [path 'h3_relval.mat']...
    }; % contains CCLE cellline names for H3 proteomics, corresponding marker ids, growth media, relative H3 proteomics
for kk = 1:numel(vars)
    load(vars{kk})
end

% impute missing values using KNN. Maybe try other imputation functions.
h3_relval = knnimpute(h3_relval);

% old variables but slightly modified
path = './../vars/';
vars = {[path 'metabolites.mat']};
for kk = 1:numel(vars)
    load(vars{kk})
end

idx = find(ismember(h3_ccle_names, celllinenames_ccle1));
tmp = length(idx);

% get relevant data
h3_relval = h3_relval(idx, :);
h3_ccle_names = h3_ccle_names(idx,1);

idx = find(ismember(celllinenames_ccle1, h3_ccle_names));
for i = 1:tmp
    model2 = model;

    % Takes in genes that are differentially expression from Z-score
    % scale
    ongenes = unique(ccleids_met(ccle_expression_metz(:,idx(i)) >= 2));
    offgenes = unique(ccleids_met(ccle_expression_metz(:,idx(i)) <= -2));
    
    % Keep the genes that match with the metabolic model.
    ongenes = intersect(ongenes, model2.rxns);
    offgenes = intersect(offgenes, model2.rxns);
    
    % set medium conditions unique to each cell line
    model2 = media(model2, h3_media(i));
    disp(i)
    
    % Get the reactions corresponding to on- and off-genes
    [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
    [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);

    % Get the flux redistribution values associated with different media component addition and deletion
    %[fluxstate_gurobi, grate_ccle_exp_dat(i,1), solverobj_ccle(i,1)] =...
    %   constrain_flux_regulation(model2, onreactions, offreactions,...
    %    kappa, rho, epsilon, mode, [], minfluxflag);

    % Add demand reactions from the metabolite list to the metabolic model
    %for m = 1:length(metabolites(:,1))
    %    tmpname = char(metabolites(m,1));
        
    % limit methionine levels for all reactions in the model; it has to be non limiting
    model3 = model2;
    [ix, pos]  = ismember({'EX_met_L(e)'}, model3.rxns);
    model3.lb(pos) = -0.5;
    rxnname = char(metabolites(:, 1)); % reaction positions of interest
    rxnpos = [find(ismember(model3.rxns, rxnname))];
    model3.c(rxnpos) = epsilon2(:,1); 

    % get the flux values from iMAT
    [fluxstate_gurobi] =  constrain_flux_regulation(model3,...
        onreactions, offreactions, kappa, rho, epsilon, mode ,[],...
        minfluxflag);
    grate_ccle_exp_dat(:,i) = fluxstate_gurobi(rxnpos);
    %model3.c(rxnpos) = 0; 
    
end

% Calculate the pearson correlation coefficients for every demand reaction
% w.r.t to H3 expression
grate_ccle_exp_dat = grate_ccle_exp_dat';
[rho, pval] = corr(grate_ccle_exp_dat, h3_relval);
rxns = metabolites(:,3);

%% Save data in Excel
filename1 = './../tables/eGEMn_prot_stats.xlsx';
colname = rxns';
rowname = h3_marks;

% Rho
xlswrite(filename1, colname, string(epsilon2), 'B1:AQ1');
xlswrite(filename1, rowname, string(epsilon2), 'A2:A22');
xlswrite(filename1, rho, string(epsilon2), 'B2:U22');

% P-value
xlswrite(filename1, colname, string(epsilon2), 'B24:AQ24');
xlswrite(filename1, rowname, string(epsilon2), 'A25:A45');
xlswrite(filename1, pval, string(epsilon2), 'B25:U45');

%% Make Figures
fig = figure;
heatmap(rho)
ax = gca;
ax.Colormap = parula;
ax.XData = h3_marks;
ax.YData = rxns;
xlabel(ax, 'Histone Markers');
ylabel(ax, 'Demand Reactions');
base = strcat('./../figures/corr/histone_mark_corr_', string(epsilon2)); 
fig_str = strcat(base, '.fig');
saveas(fig, fig_str);
end 