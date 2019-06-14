%% @author: Scott Campit
function [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = metabolic_sensitivity(model, compartment,...
    epsilon2, scaling, exp)
% metabolic_sensitivity.m displays the values corresponding to several demand
% reactions and excess/depletion of a specific medium component.

% INPUTS:
    % model: genome scale metabolic model
    % rxnpos: list of reactions that will be used in the analysis
    % epsilon2: a weight that will be used for the secondary obcomponentective coefficient
    % scaling: if the data should be scaled/normalized

% OUTPUTS:
    % excess_flux: metabolic flux corresponding to excess medium
    % depletion_flux: ''' depletion medium
    % excess_redcost: reduced cost corresponding to excess medium
    % depletion_redcost: ''' depletion medium
    % excess_shadow: shadow costs corresponding to excess medium
    % depletion_shadow: ''' depletion medium

%% metabolic_sensitivity.m
load ./../vars/supplementary_software_code media_exchange1
var = {'./../vars/metabolites.mat', './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk});
end

% init some params
minfluxflag = 0;
posgluc = 1385;  % glucose uptake reaction in RECON1
biomassobcomponentpos = 3743; % biomass rxn position
rxnname = char(metabolites(:, 1)); % reaction positions of interest

% Switch between modules for different types of analyses
switch exp
    % Case 1: Use the most dynamic range of metabolic fluxes
    case 'fba'  
        % Kappa values are weights for both excess and depletion of medium.
        % kappa = 10 is excess, and kappa = 0.01 is depletion
        for kappatype = 1:2
            if kappatype == 1
                kappa  = 10;
            elseif kappatype == 2
                kappa = 0.01;
            end
            
            % Each medium component has different uptake rates
            for component = 1:length(mediareactions1(:, 1)) % 50 times
                if (kappatype == 2) & (ismember(component,[2,3,5:19])) % trace elements
                    weight = kappa/100;
                elseif (kappatype == 1) & (ismember(component,[1;4])) % glucose or glutamine
                    weight = 3;
                end
                
                model2 = model;

                % Make the methylation exchange reaction have a fixed LB of
                % -0.5 to be non-limiting
                [~, pos] = ismember({'EX_met_L(e)'}, model2.rxns);
                model2.lb(pos) = -0.5;
                
                % Get the reaction position for the various medium components
                % and constrain the lower bound with a value proportional to
                % the weight
                [~, pos]  = ismember(mediareactions1(component,1), model2.rxns);
                model2.lb(pos) = -1*media_exchange1(component,1)*weight;
                %model2 = media(model2, 'RPMI');
                %model2.c(biomassobcomponentpos) = 1;

                %% solve using FBA (cobratoolbox) to get growth rate
                soln = optimizeCbModel(model2);

                % Metabolic fluxes for a single reaction
                grate_str = ['media_xchange_grate_', num2str(kappatype),...
                    '(component,1:20) = soln.v(biomassobcomponentpos);'];

                % Reduced costs for a single reaction
                rc_str = ['media_xchange_grate_rc_', num2str(kappatype),...
                    '(component,1:20) = soln.w(biomassobcomponentpos);'];

                if ~isempty(soln.v) & ~isnan(soln.v)
                    eval(grate_str)
                end
                if ~isempty(soln.w) & ~isnan(soln.w)
                    eval(rc_str)
                end

                if kappatype == 1
                    excess_xgrate(component) = media_xchange_grate_1(component);
                    excess_xredcost(component) = media_xchange_grate_rc_1(component);
                   
                end
                if kappatype == 2
                    depletion_xgrate(component) = media_xchange_grate_2(component);
                    depletion_xredcost(component) = media_xchange_grate_rc_2(component);
                    
                end
                model2.c(biomassobcomponentpos) = 0;
                
                %% Obtain flux values when using epsilon2 as the obcomponentective coefficient for the reaction of interest.
                model3 = model2;
                rxnpos = [find(ismember(model3.rxns, rxnname))];
                
                % epsilon values for excess and depletion respectively
                if kappatype == 1
                    model3.c(rxnpos) = epsilon2(:, 1);
                elseif kappatype == 2
                    model3.c(rxnpos) = epsilon2(:, 2);
                end
                
                soln = optimizeCbModel(model3);

                % Metabolic fluxes
                flux_str = ['media_xchange_rxn_flux_', num2str(kappatype),...
                    '(component,1:20) = soln.v(rxnpos);'];
                % Reduced costs
                rc_str = ['media_xchange_rxn_rc_', num2str(kappatype),...
                    '(component,1:20) = soln.w(rxnpos);'];
                % Shadow prices
                for met = 1:length(metabolites(:,2))
                    tmp_met{met} = [char(metabolites(met, 2)) '[' compartment ']'];
                end
                tmp_met = tmp_met';
                
                met_pos = find(ismember(model2.mets, tmp_met));
                sp_str = ['media_xchange_rxn_sp_', num2str(kappatype),...
                    '(component,:) = soln.y(met_pos);'];

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
                    excess_flux = media_xchange_rxn_flux_1;
                    excess_redcost = media_xchange_rxn_rc_1;
                    excess_shadow = media_xchange_rxn_sp_1;
                end
                if kappatype == 2
                    depletion_flux = media_xchange_rxn_flux_2;
                    depletion_redcost = media_xchange_rxn_rc_2;
                    depletion_shadow = media_xchange_rxn_sp_2;
                end
                model3.c(rxnpos) = 0;
                disp(component)
            end
            disp(kappatype)
        end
    
    % Run the code using the FVA scheme for histone markers only
    case 'fva'
        for kappatype = 1:2
            if kappatype == 1
                kappa  = 10;
            elseif kappatype == 2
                kappa = 0.01;
            end
            
            % Each medium component has different uptake rates
            for component = 1:length(mediareactions1(:, 1)) % 50 times
                if (kappatype == 2) & (ismember(component,[2,3,5:19])) % trace elements
                    weight = kappa/100;
                elseif (kappatype == 1) & (ismember(component,[1;4])) % glucose or glutamine
                    weight = 3;
                end
                
                model2 = model;

                % Make the methylation exchange reaction have a fixed LB of
                % -0.5 to be non-limiting
                [~, pos] = ismember({'EX_met_L(e)'}, model2.rxns);
                model2.lb(pos) = -0.5;

                % Get the reaction position for the various medium components
                % and constrain the lower bound with a value proportional to
                % the weight
                [~, pos]  = ismember(mediareactions1(component,1), model2.rxns);
                model2.lb(pos) = -1*media_exchange1(component,1)*weight;
                model2.c(biomassobcomponentpos) = 1;

                %% solve using FBA (cobratoolbox) to get growth rate
                soln = optimizeCbModel(model2);

                % Metabolic fluxes for a single reaction
                grate_str = ['media_xchange_grate_', num2str(kappatype),...
                    '(component,1:20) = soln.v(biomassobcomponentpos);'];

                % Reduced costs for a single reaction
                rc_str = ['media_xchange_grate_rc_', num2str(kappatype),...
                    '(component,1:20) = soln.w(biomassobcomponentpos);'];

                if ~isempty(soln.v) & ~isnan(soln.v)
                    eval(grate_str)
                end
                if ~isempty(soln.w) & ~isnan(soln.w)
                    eval(rc_str)
                end

                if kappatype == 1
                    excess_xgrate(component) = media_xchange_grate_1(component);
                    excess_xredcost(component) = media_xchange_grate_rc_1(component);
                   
                end
                if kappatype == 2
                    depletion_xgrate(component) = media_xchange_grate_2(component);
                    depletion_xredcost(component) = media_xchange_grate_rc_2(component);
                    
                end
                model2.c(biomassobcomponentpos) = 0;
                
                %% Obtain flux values when using epsilon2 as the obcomponentective coefficient for the reaction of interest.
                model3 = model2;
                rxnpos = [find(ismember(model3.rxns, rxnname))];
                rxn_nam = model3.rxns(rxnpos);
                
                % get histone marks
                histmark = model3.rxns([3786:3789]);
                
                % epsilon values for excess and depletion respectively
                if kappatype == 1
                    model3.c(rxnpos) = epsilon2(:, 1);
                elseif kappatype == 2
                    model3.c(rxnpos) = epsilon2(:, 2);
                end
                
                % Run FVA for methylation and acetylation reactions.
                [minFlux, maxFlux] = fluxVariability(model3, 100, 'max', rxn_nam);

                % Metabolic fluxes
                flux_str = ['media_xchange_rxn_flux_', num2str(kappatype),...
                    '(component,1:20) = maxFlux;'];

                if ~isempty(soln.v) & ~isnan(soln.v)
                    eval(flux_str)
                end
          
                if kappatype == 1
                    excess_flux = media_xchange_rxn_flux_1;
                end
                if kappatype == 2
                    depletion_flux = media_xchange_rxn_flux_2;
                end
                model3.c(rxnpos) = 0;
                disp(component)
            end
            disp(kappatype)
        end
end

%% Scaling the data for visualization purposes
excess_flux = zscore(excess_flux);
excess_redcost = zscore(excess_redcost);
excess_shadow = zscore(excess_shadow);

depletion_flux = zscore(depletion_flux);
depletion_redcost = zscore(depletion_redcost);
depletion_shadow = zscore(depletion_shadow);

%% Replace 0 with NaN
excess_flux(excess_flux==0) = NaN;
excess_redcost(excess_redcost==0) = NaN;
excess_shadow(excess_shadow==0) = NaN;

depletion_flux(depletion_flux==0) = NaN;
depletion_redcost(depletion_redcost==0) = NaN;
depletion_shadow(depletion_shadow==0) = NaN;

%% Save metabolic flux data as excel file
% input filename for saving
filename1 = './../tables/eGEMn.xlsx';
colname = metabolites(:, 3)';
rowname = mediareactions1(:, 2);

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

%% Barplot figures
% Prepare figure labels and variables
medium_labels = mediareactions1(:,2);
reaction_labels = metabolites(:,3);

fig1 = figure;
subplot(1,2,1);
heatmap(excess_flux)
ax1 = gca;
ax1.XData = reaction_labels;
ax1.YData = medium_labels;
ax1.Title = 'Metabolic flux in excess medium';
xlabel(ax1, 'Demand reactions');
ylabel(ax1, 'Medium component');

subplot(1,2,2);
heatmap(depletion_flux)
ax2 = gca;
ax2.XData = reaction_labels;
ax2.YData = medium_labels;
ax2.Title = 'Metabolic flux in depleted medium';
xlabel(ax2, 'Demand reactions');
ylabel(ax2, 'Medium component');



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
subplot(2,3,1);
heatmap(excess_flux)
ax1 = gca;
ax1.XData = reaction_labels;
ax1.YData = medium_labels;
ax1.Title = 'Metabolic flux in excess medium';
xlabel(ax1, 'Demand reactions');
ylabel(ax1, 'Medium component');

subplot(2,3,4);
heatmap(depletion_flux)
ax2 = gca;
ax2.XData = reaction_labels;
ax2.YData = medium_labels;
ax2.Title = 'Metabolic flux in depleted medium';
xlabel(ax2, 'Demand reactions');
ylabel(ax2, 'Medium component');

subplot(2,3,2);
heatmap(excess_shadow)
ax3 = gca;
ax3.XData = reaction_labels;
ax3.YData = medium_labels;
ax3.Title = 'Shadow price in excess medium';
xlabel(ax3, 'Demand reactions');
ylabel(ax3, 'Medium component');

subplot(2,3,5);
heatmap(depletion_shadow)
ax4 = gca;
ax4.XData = reaction_labels;
ax4.YData = medium_labels;
ax4.Title = 'Shadow price in depleted medium';
xlabel(ax4, 'Demand reactions');
ylabel(ax4, 'Medium component');

subplot(2,3,3);
heatmap(excess_redcost)
ax5 = gca;
ax5.XData = reaction_labels;
ax5.YData = medium_labels;
ax5.Title = 'Reduced cost in excess medium';
xlabel(ax5, 'Demand reactions');
ylabel(ax5, 'Medium component');

subplot(2,3,6);
heatmap(depletion_redcost)
ax6 = gca;
ax6.XData = reaction_labels;
ax6.YData = medium_labels;
ax6.Title = 'Reduced cost in depleted medium';
xlabel(ax6, 'Demand reactions');
ylabel(ax6, 'Medium component');

base = strcat('./../figures/new-model/eGEMn_', string(epsilon2));
fig1_str = strcat(base, '.fig');

saveas(fig1, fig1_str);

% Create a heatmap for the growth rates in excess and depleted medium
fig2 = figure;

dat = [excess_xgrate; depletion_xgrate];
dat = dat';

heatmap(dat)
ax1 = gca;
%ax1.XData = ['Excess Medium'; 'Depleted Medium'];
ax1.YData = medium_labels;
ax1.Title = 'Growth rate in excess medium';
xlabel(ax1, 'Growth Rates');
ylabel(ax1, 'Medium component');

subplot(1,2,1);
heatmap(excess_xgrate)
ax1 = gca;
ax1.XData = 'Growth Rate';
ax1.YData = medium_labels;
ax1.Title = 'Growth rate in excess medium';
xlabel(ax1, 'Demand reactions');
ylabel(ax1, 'Medium component');

subplot(1,2,2);
heatmap(depletion_xgrate)
ax2 = gca;
ax2.XData = 'Growth Rate';
ax2.YData = medium_labels;
ax2.Title = 'Growth rate in depleted medium';
xlabel(ax2, 'Demand reactions');
ylabel(ax2, 'Medium component');

base = strcat('./../figures/new-model/eGEMn_grate_', string(epsilon2));
fig2_str = strcat(base, '.fig');
saveas(fig2, fig2_str);

%clear
end
