%% histone_corr calculates the correlation value between various histone
...markers and the metabolic flux obtained from the iMAT algorithm
function [grate_ccle_exp_dat, rho, pval] = histone_corr(model, dat, cell_names, reactions_of_interest,...
compartment, medium, epsilon2, mode, epsilon, rho, kappa, minfluxflag, type,...
fva_grate)

%% INPUTS:
    % model: A structure representing the initial genome scale model
    % dat: A numerical array containing proteomics data used for the correlation
    % reactions_of_interest: A cell array containing reactions that will be studied
    % compartment: A character specifying the subcellular compartment of interest based on BiGG classifications
    % mode: describes the input list for constrain flux regulation
        % 0 =  genes
        % 1 = reactions
    % epsilon: The parameter for the minimum flux of on-reactions for constrain flux regulation. The default value is 1E-3
    % epsilon2: The array of objective coefficients for reaction of interest
    % rho: The relative weight for on-reactions. Used as a parameter for constrain flux regulation
    % kappa: The relative weight for off-reactions. Used as a parameter for constrain flux regulation
    % minfluxflag: Binary input for parsimonious flux balance analysis
        % 0 = no pFBA
        % 1 = pFBA

%% OUTPUTS:
    % correl: A numerical array of correlation values associated with each histone marker/rxn
    % pval: A numerical array of the p-value associated from the pearson correlation computation
    % cell_line_match: A cell array containing the cell lines that matched between gene expression and proteomics

%% histone_corr
load ./../vars/supplementary_software_code...
    celllinenames_ccle1... % CCLE cellline names
    ccleids_met... % Gene symbols
    ccle_expression_metz % Z-transformed gene expression

BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); % biomass rxn position in eGEM model
obj_coefs = epsilon2{:};
%obj_coefs = cell2mat(obj_coefs);

% impute missing values using KNN and scale from [0,1]
dat = knnimpute(dat);
dat = normalize(dat, 'range');

% Match data from gene expression and histone proteomics to get proteomics
% data that will be used downstream
idx = find(ismember(cell_names, celllinenames_ccle1));
tmp = length(idx);
dat = dat(idx, :);
cell_names = cell_names(idx,1);

% Change idx to map to gene expression array and iterate for all 885 cancer
% cell lines that match between genexp and proteomics dataset
idx = find(ismember(celllinenames_ccle1, cell_names));
for i = 1:tmp
    disp(i)
    model2 = model;

    ongenes = unique(ccleids_met(ccle_expression_metz(:,idx(i)) >= 2));
    offgenes = unique(ccleids_met(ccle_expression_metz(:,idx(i)) <= -2));
    ongenes = intersect(ongenes, model2.rxns);
    offgenes = intersect(offgenes, model2.rxns);
    [ix, pos]  = ismember({'EX_met_L(e)'}, model2.rxns);
    model2 = media(model2, medium);
    model2.lb(pos) = -0.5;
    model2.c(BIOMASS_OBJ_POS) = 1;

    [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
    [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);

    % Get the demand reaction positions of interest and calculate metabolic
    % flux for each cell line using the iMAT algorithm
    rxnname = char(reactions_of_interest(:, 1));
    switch type
        case 'non-competitive_cfr'
            for rxn = 1:length(reactions_of_interest)
                disp(rxn)
                rxnpos = [find(ismember(model2.rxns, reactions_of_interest(rxn)))];
                model2.c(rxnpos) = obj_coefs(rxn, 1);
                [flux, ~, ~] = constrain_flux_regulation(model2,  ...
                    onreactions, offreactions, kappa, rho, epsilon, mode, [], ...
                    minfluxflag);
                grate_ccle_exp_dat(i,rxn) = flux(rxnpos);
            end
        case 'competitive_cfr'
            rxnpos = [find(ismember(model2.rxns, rxnname))];
            model2.c(rxnpos) = obj_coefs(:, 1);
            [flux, ~, ~] =  constrain_flux_regulation(model2,...
                onreactions, offreactions, kappa, rho, epsilon, mode , [], ...
                minfluxflag);
            grate_ccle_exp_dat(i,:) = flux(rxnpos);
        case 'fva'
            rxnpos = [find(ismember(model2.rxns, rxnname))];
            model2.c(rxnpos) = obj_coefs(:, 1);
            [~, ~, ~, ~, flux, ~] =...
                calc_metabolic_metrics(model2, rxnpos, [], fva_grate,...
                'max', reactions_of_interest, obj_coefs(:, 1), type);
            grate_ccle_exp_dat(i,:) = flux(rxnpos);
end

% Calculate the pearson correlation coefficients for every demand reaction
disp(size(grate_ccle_exp_dat))
disp(size(dat))
[rho, pval] = corr(grate_ccle_exp_dat, dat);
rxns = metabolites(:,3);

%% Save data in Excel
filename1 = './../tables/eGEMn_allDM.xlsx';
colname = rxns';
rowname = h3_marks;
xlswrite(filename1, colname, string(obj_coefs), 'B1:U1');
xlswrite(filename1, rowname, string(obj_coefs), 'A2:A51');
xlswrite(filename1, excess_flux, string(obj_coefs), 'B2:U51');

%% Make Figures
fig = figure;
heatmap(rho)
ax = gca;
ax.Colormap = parula;
ax.XData = h3_marks;
ax.YData = rxns;
ax.Title = 'Histone markers and metabolic flux correlation'
xlabel(ax, 'Histone Markers');
ylabel(ax, 'Cancer Cell Lines (CCLE)');
saveas(fig, ['./../figures/corr/histone_mark_corr.fig']);
end
