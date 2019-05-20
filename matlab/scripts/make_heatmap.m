%% @author: Scott Campit
function [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = make_heatmap(model, compartment,...
    epsilon2, scaling)
% make_heatmap.m displays the values corresponding to several demand 
% reactions and excess/depletion of a specific medium component.

% INPUTS:
    % model: genome scale metabolic model
    % rxnpos: list of reactions that will be used in the analysis
    % epsilon2: a weight that will be used for the secondary objective coefficient 
    % scaling: if the data should be scaled/normalized 

% OUTPUTS:
    % excess_flux: metabolic flux corresponding to excess medium
    % depletion_flux: ''' depletion medium
    % excess_redcost: reduced cost corresponding to excess medium
    % depletion_redcost: ''' depletion medium
    % excess_shadow: shadow costs corresponding to excess medium
    % depletion_shadow: ''' depletion medium
    
%% make_heatmap.m
load supplementary_software_code media_exchange1 
var = {'./../vars/metabolites.mat', './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk})
end

minfluxflag = 0;
posgluc = 1385;  % glucose uptake reaction in RECON1
biomassobjpos = find(model.c); % biomass objective

% Add demand reactions from the metabolite list to the metabolic model
for m = 1:length(metabolites(:,1))
    tmp_met = char(metabolites(m,2));
    tmp = [tmp_met '[' compartment '] -> '];
    tmpname = char(metabolites(m,1));
    model2 = addReaction(model, tmpname, 'reactionFormula', tmp);
    
    % Kappa values are weights for both excess and depletion of medium.
    % kappa = 10 is excess, and kappa = 0.01 is depletion
    for kappatype = 1:2
        if kappatype == 1 
            kappa  = 10; 
        elseif kappatype == 2
            kappa = 0.01;
        end

        for j = 1:length(mediareactions1(:, 1))
            kappa1 = kappa;
            if (kappatype == 2) & (ismember(j,[2,3,5:19])) % trace elements
                kappa1 = kappa/100;
            elseif (kappatype == 1) & (ismember(j,[1;4])) % glucose or glutamine
                kappa1 = 3;
            end

            % Make the methylation exchange reaction have a fixed LB of
            % -0.5 to be non-limiting
            [~, pos] = ismember({'EX_met_L(e)'}, model2.rxns);
            model2.lb(pos) = -0.5;    

            % Get the reaction position for the various medium components
            % and constrain the lower bound with a value proportional to
            % the weight kappa1
            [~, pos]  = ismember(mediareactions1(j,1), model2.rxns);
            model2.lb(pos) = -media_exchange1(j,1)*kappa1;
            
            %% solve using FBA (cobratoolbox) to get growth rate
            soln = optimizeCbModel(model2);
            
            % Metabolic fluxes for a single reaction
            grate_str = ['media_xchange_grate_', num2str(kappatype),...
                '(j,1) = soln.v(biomassobjpos);'];
            
            % Reduced costs for a single reaction
            rc_str = ['media_xchange_rc_', num2str(kappatype),...
                '(j,1) = soln.w(biomassobjpos);'];
            
            % Shadow prices for a single metabolite
            tmp_met = [char(metabolites(m,2)) '[' compartment ']'];
            met_pos = find(ismember(model2.mets,tmp_met));
            sp_str = ['media_xchange_grate_', num2str(kappatype),...
                '(j,1) = soln.y(met_pos);'];
            
            if ~isempty(soln.v) & ~isnan(soln.v)
                eval(grate_str)
            end
            if ~isempty(soln.w) & ~isnan(soln.w)
                eval(rc_str)
            end
            if ~isempty(soln.y) & ~isnan(soln.y)
                eval(sp_str)
            end
            
            if kappatype == 1
                excess_xgrate(j,m) = media_xchange_grate_1(j,1);
                excess_xredcost(j,m) = media_xchange_rc_1(j,1);
                excess_xshadow(j,m) = media_xchange_sp_1(j,1);
            end
            if kappatype == 2
                depletion_xgrate(j,m) = media_xchange_grate_2(j,1);
                depletion_xredcost(j,m) = media_xchange_rc_2(j,1);
                depletion_xshadow(j,m) = media_xchange_sp_2(j,1);
            end
            
            %% Obtain flux values when using epsilon2 as the objective coefficient for the reaction of interest.
            model3 = model2;
            rxnpos = [find(ismember(model3.rxns,tmpname))];
            model3.c(rxnpos) = epsilon2;
            soln = optimizeCbModel(model3);

            % Metabolic fluxes
            flux_str = ['media_xchange_rxn_flux_', num2str(kappatype),...
                '(j,1) = soln.v(rxnpos);'];
            % Reduced costs
            rc_str = ['media_xchange_rxn_rc_', num2str(kappatype),...
                '(j,1) = soln.w(rxnpos);'];
            % Shadow prices
            tmp_met = [char(metabolites(m,2)) '[' compartment ']'];
            met_pos = find(ismember(model3.mets,tmp_met));
            sp_str = ['media_xchange_rxn_sp_', num2str(kappatype),...
                '(j,1) = soln.y(met_pos);'];

            if ~isempty(soln.v) & ~isnan(soln.v)
                eval(flux_str)
            end
            if ~isempty(soln.w) & ~isnan(soln.w)
                eval(rc_str)
            end
            if ~isempty(soln.y) & ~isnan(soln.y)
                eval(sp_str)
            end
            
            if kappatype == 1
                excess_flux(j,m) = media_xchange_rxn_flux_1(j,1);
                excess_redcost(j,m) = media_xchange_rxn_rc_1(j,1);
                excess_shadow(j,m) = media_xchange_rxn_sp_1(j,1);
            end
            if kappatype == 2
                depletion_flux(j,m) = media_xchange_rxn_flux_2(j,1);
                depletion_redcost(j,m) = media_xchange_rxn_rc_2(j,1);
                depletion_shadow(j,m) = media_xchange_rxn_sp_2(j,1);
            end
            disp(j)
        end
        disp(kappatype)
    end
end

%% Get flux difference
%exces_flux = exces_flux - excess_xchange_flux;
%depletion_flux = depletion_flux - depletion_xchange_flux;
%excess_shadow = excess_shadow - excess_xchange_shadow;
%depletion_shadow = depletion_shadow - depletion_shadow;
%excess_redcost = excess_redcost - excess_redcost;
%depletion_redcost = depletion_redcost - depletion_redcost;

%% Save metabolic flux data as excel file
% input filename for saving
filename1 = './../tables/eGEMn.xlsx';
colname = metabolites(:,3)';
rowname = mediareactions1(:,2);

% Excess flux
xlswrite(filename1, colname, string(epsilon2), 'B1:U1');
xlswrite(filename1, rowname, string(epsilon2), 'A2:A51');
xlswrite(filename1, excess_flux, string(epsilon2), 'B2:U51');
% Depleted flux
xlswrite(filename1, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename1, rowname, string(epsilon2), 'A2:A51');
xlswrite(filename1, depletion_flux, string(epsilon2), 'X2:AQ51');
% Excess reduced cost
xlswrite(filename1, colname, string(epsilon2), 'B1:U1');
xlswrite(filename1, rowname, string(epsilon2), 'A54:A103');
xlswrite(filename1, excess_redcost, string(epsilon2), 'B54:U103');
% Depleted reduced cost
xlswrite(filename1, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename1, rowname, string(epsilon2), 'A54:A103');
xlswrite(filename1, depletion_redcost, string(epsilon2), 'X54:AQ103');
% Excess shadow price
xlswrite(filename1, colname, string(epsilon2), 'B1:U1');
xlswrite(filename1, rowname, string(epsilon2), 'A106:A155');
xlswrite(filename1, excess_shadow, string(epsilon2), 'B106:U155');
% Depleted shadow price
xlswrite(filename1, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename1, rowname, string(epsilon2), 'A106:A155');
xlswrite(filename1, depletion_shadow, string(epsilon2), 'X106:AQ155');

%% Save grate data as excel file
% input filename for saving
filename2 = './../tables/eGEMn_grate.xlsx';
colname = metabolites(:,3)';
rowname = mediareactions1(:,2);

% Excess grate
xlswrite(filename2, colname, string(epsilon2), 'B1:U1');
xlswrite(filename2, rowname, string(epsilon2), 'A2:A51');
xlswrite(filename2, excess_xgrate, string(epsilon2), 'B2:U51');
% Depleted grate
xlswrite(filename2, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename2, rowname, string(epsilon2), 'A2:A51');
xlswrite(filename2, depletion_xgrate, string(epsilon2), 'X2:AQ51');
% Excess reduced cost
xlswrite(filename2, colname, string(epsilon2), 'B1:U1');
xlswrite(filename2, rowname, string(epsilon2), 'A54:A103');
xlswrite(filename2, excess_xredcost, string(epsilon2), 'B54:U103');
% Depleted reduced cost
xlswrite(filename2, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename2, rowname, string(epsilon2), 'A54:A103');
xlswrite(filename2, depletion_xredcost, string(epsilon2), 'X54:AQ103');
% Excess shadow price
xlswrite(filename2, colname, string(epsilon2), 'B1:U1');
xlswrite(filename2, rowname, string(epsilon2), 'A106:A155');
xlswrite(filename2, excess_xshadow, string(epsilon2), 'B106:U155');
% Depleted shadow price
xlswrite(filename2, colname, string(epsilon2), 'X1:AQ1');
xlswrite(filename2, rowname, string(epsilon2), 'A106:A155');
xlswrite(filename2, depletion_xshadow, string(epsilon2), 'X106:AQ155');

%% Turn 0 -> NaN for heatmap visualization
excess_flux(excess_flux == 0) = NaN;
depletion_flux(depletion_flux == 0) = NaN;
excess_shadow(excess_shadow == 0) = NaN;
depletion_shadow(depletion_shadow == 0) = NaN;
excess_redcost(excess_redcost == 0) = NaN;
depletion_redcost(depletion_redcost == 0) = NaN;

%% Heatmap figures
% Still needs work:
    % Make `excess` green and `depletion` red
    % Color gradient: metabolic flux > reduced costs > shadow prices
    % Incorporate plotly 
    % Is there a way to dynamically size .pngs?
    
% Prepare figure labels and variables
medium_labels = mediareactions1(:,2);
reaction_labels = metabolites(:,3);

fig1 = figure;
vars = [excess_flux, depletion_flux, excess_shadow, depletion_shadow,...
    excess_redcost, depletion_redcost];
title_label = ['Metabolic flux in excess medium',...
    'Metabolic flux in depleted medium',...
    'Shadow price in excess medium',...
    'Shadow price in depleted medium',...
    'Reduced cost in excess medium',...
    'Reduced cost in depleted medium'
    ];
for i = 1:length(vars)
    subplot(2,3,i);
    heatmap(vars(i))
    ax = ['ax' i];
    ax = gca;
    ax.XData = reaction_labels;
    ax.YData = medium_labels;
    ax.Title = title_label(i);
    xlabel(ax, 'Demand reactions');
    ylabel(ax, 'Medium component');

base = strcat('./../figures/new-model/eGEMn_', string(epsilon2));
fig1_str = strcat(base, '.fig');
saveas(fig1, fig1_str);

% Create a heatmap for the growth rates in excess and depleted medium
fig2 = figure;
subplot(1,2,1);
heatmap(excess_xgrate)
ax1 = gca;
ax1.XData = reaction_labels;
ax1.YData = medium_labels;
ax1.Title = 'Growth rate in excess medium';
xlabel(ax1, 'Demand reactions');
ylabel(ax1, 'Medium component');

subplot(1,2,2);
heatmap(depletion_xgrate)
ax2 = gca;
ax2.XData = reaction_labels;
ax2.YData = medium_labels;
ax2.Title = 'Growth rate in depleted medium';
xlabel(ax2, 'Demand reactions');
ylabel(ax2, 'Medium component');

base = strcat('./../figures/new-model/eGEMn_grate_', string(epsilon2));
fig2_str = strcat(base, '.fig');
saveas(fig2, fig2_str);

end