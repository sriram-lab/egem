%% make_heatmap.m
function [excess_flux, depletion_flux, excess_redcost, depletion_redcost,...
    excess_shadow, depletion_shadow] = make_heatmap(model, compartment,...
    epsilon2, scaling)
% make_heatmap.m displays an array of correlation values corresponding to
% several demand reactions and the quantity of the histone maker based on
% metabolic perturbation. 

% INPUTS:
    % model: genome scale metabolic model
    % rxnpos: list of reactions that will be used in the analysis
    % epsilon2: a weight that will be used for the secondary objective coefficient 
    % scaling: if the data should be scaled/normalized 
    
%% 
load supplementary_software_code media_exchange1 
var = {'metabolites.mat', 'cellmedia.mat', 'mediareactions1.mat'};
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
    %tmpname = 'EX_KAC';
    model2 = addReaction(model, tmpname, 'reactionFormula', tmp);
    
    % Kappa values are weights for both excess and depletion of medium.
    % kappa = 10 is excess, and kappa = 0.01 is depletion
    for kappatype = 1:2
        if kappatype == 1 
            kappa  = 10; 
        elseif kappatype == 2
            kappa = 0.01;
        end

        for j = 1:length(mediareactions1(:,1))
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
            
            %% solve using FBA (cobratoolbox) without setting model.c(rxn) = epsilon2 
            soln = optimizeCbModel(model2);
            
            % Metabolic fluxes for a single reaction
            flux_str = ['media_xchange_flux_', num2str(kappatype),...
                '(j,1) = soln.v(biomassobjpos);'];
            
            % Reduced costs for a single reaction
            rc_str = ['media_xchange_rc_', num2str(kappatype),...
                '(j,1) = soln.w(biomassobjpos);'];
            
            % Shadow prices for a single metabolite
            tmp_met = [char(metabolites(m,2)) '[' compartment ']'];
            met_pos = find(ismember(model2.mets,tmp_met));
            sp_str = ['media_xchange_sp_', num2str(kappatype),...
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
                excess_xchange_flux(j,m) = media_xchange_flux_1(j,1);
                excess_xchange_redcost(j,m) = media_xchange_rc_1(j,1);
                excess_xchange_shadow(j,m) = media_xchange_sp_1(j,1);
            end
            if kappatype == 2
                depletion_xchange_flux(j,m) = media_xchange_flux_2(j,1);
                depletion_xchange_redcost(j,m) = media_xchange_rc_2(j,1);
                depletion_xchange_shadow(j,m) = media_xchange_sp_2(j,1);
            end
            
            %% Obtain flux values when using epsilon2 as the objective
            % coefficient for the reaction of interest.
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
                excess_xchange_rxn_flux(j,m) = media_xchange_rxn_flux_1(j,1);
                excess_xchange_rxn_redcost(j,m) = media_xchange_rxn_rc_1(j,1);
                excess_xchange_rxn_shadow(j,m) = media_xchange_rxn_sp_1(j,1);
            end
            if kappatype == 2
                depletion_xchange_rxn_flux(j,m) = media_xchange_rxn_flux_2(j,1);
                depletion_xchange_rxn_redcost(j,m) = media_xchange_rxn_rc_2(j,1);
                depletion_xchange_rxn_shadow(j,m) = media_xchange_rxn_sp_2(j,1);
            end
            disp(j)
        end
        disp(kappatype)
    end
end

%% Normalize the flux, shadow price, and reduced costs.
excess_xchange_rxn_flux = excess_xchange_rxn_flux - excess_xchange_flux;
depletion_xchange_rxn_flux = depletion_xchange_rxn_flux - depletion_xchange_flux;
excess_xchange_rxn_shadow = excess_xchange_rxn_shadow - excess_xchange_shadow;
depletion_xchange_rxn_shadow = depletion_xchange_rxn_shadow - depletion_xchange_rxn_shadow;
excess_xchange_rxn_redcost = excess_xchange_rxn_redcost - excess_xchange_rxn_redcost;
depletion_xchange_rxn_redcost = depletion_xchange_rxn_redcost - depletion_xchange_rxn_redcost;

%% Turn 0 -> NaN for heatmap visualization
excess_xchange_rxn_flux(excess_xchange_rxn_flux == 0) = NaN;
depletion_xchange_rxn_flux(depletion_xchange_rxn_flux == 0) = NaN;
excess_xchange_rxn_shadow(excess_xchange_rxn_shadow == 0) = NaN;
depletion_xchange_rxn_shadow(depletion_xchange_rxn_shadow == 0) = NaN;
excess_xchange_rxn_redcost(excess_xchange_rxn_redcost == 0) = NaN;
depletion_xchange_rxn_redcost(depletion_xchange_rxn_redcost == 0) = NaN;

%% Scaling arguments <!-- Only takes into account the model with epsilon2 -->
if (~exist('scaling', 'var')) || (isempty(scaling))
    excess_flux = excess_xchange_rxn_flux; 
    depletion_flux = depletion_xchange_rxn_flux;
    excess_shadow = excess_xchange_rxn_shadow;
    depletion_shadow = depletion_xchange_rxn_shadow;
    excess_redcost = excess_xchange_rxn_redcost;
    depletion_redcost = depletion_xchange_rxn_redcost;
elseif scaling == 'minmax'
    excess_flux = normalize(excess_xchange_rxn_flux, 'range');
    depletion_flux = normalize(depletion_xchange_rxn_flux, 'range');
    excess_shadow = normalize(excess_xchange_rxn_shadow, 'range');
    depletion_shadow = normalize(depletion_xchange_rxn_shadow, 'range');
    excess_redcost = normalize(excess_xchange_rxn_redcost, 'range');
    depletion_redcost = normalize(depletion_xchange_rxn_redcost, 'range');
elseif scaling == 'sttdev'
    excess_flux = normalize(excess_xchange_rxn_flux, 'scale');
    depletion_flux = normalize(depletion_xchange_rxn_flux, 'scale');
    excess_shadow = normalize(excess_xchange_rxn_shadow, 'scale');
    depletion_shadow = normalize(depletion_xchange_rxn_shadow, 'scale');
    excess_redcost = normalize(excess_xchange_rxn_redcost, 'scale');
    depletion_redcost = normalize(depletion_xchange_rxn_redcost, 'scale');
elseif scaling == '1-norm'
    excess_flux = normalize(excess_xchange_rxn_flux, 'norm', 1);
    depletion_flux = normalize(depletion_xchange_rxn_flux, 'norm', 1);
    excess_shadow = normalize(excess_xchange_rxn_shadow, 'norm', 1);
    depletion_shadow = normalize(depletion_xchange_rxn_shadow, 'norm', 1);
    excess_redcost = normalize(excess_xchange_rxn_redcost, 'norm', 1);
    depletion_redcost = normalize(depletion_xchange_rxn_redcost, 'norm', 1);
elseif scaling == '2-norm'
    excess_flux = normalize(excess_xchange_rxn_flux, 'norm', 2);
    depletion_flux = normalize(depletion_xchange_rxn_flux, 'norm', 2);
    excess_shadow = normalize(excess_xchange_rxn_shadow, 'norm', 2);
    depletion_shadow = normalize(depletion_xchange_rxn_shadow, 'norm', 2);
    excess_redcost = normalize(excess_xchange_rxn_redcost, 'norm', 2);
    depletion_redcost = normalize(depletion_xchange_rxn_redcost, 'norm', 2);
end

%% Heatmap figures

% Prepare figure labels and variables
medium_labels = mediareactions1(:,2);
reaction_labels = metabolites(:,3);

fig = figure;

subplot(3,2,1);
heatmap(excess_flux)
ax1 = gca;
ax1.XData = reaction_labels;
ax1.YData = medium_labels;
ax1.Title = 'Metabolic flux in excess medium';
xlabel(ax1, 'Demand reactions');
ylabel(ax1, 'Medium component');

subplot(3,2,2);
heatmap(depletion_flux)
ax2 = gca;
ax2.XData = reaction_labels;
ax2.YData = medium_labels;
ax2.Title = 'Metabolic flux in depleted medium';
xlabel(ax2, 'Demand reactions');
ylabel(ax2, 'Medium component');

subplot(3,2,3);
heatmap(excess_shadow)
ax3 = gca;
ax3.XData = reaction_labels;
ax3.YData = medium_labels;
ax3.Title = 'Shadow price in excess medium';
xlabel(ax3, 'Demand reactions');
ylabel(ax3, 'Medium component');

subplot(3,2,4);
heatmap(depletion_shadow)
ax4 = gca;
ax4.XData = reaction_labels;
ax4.YData = medium_labels;
ax4.Title = 'Shadow price in depleted medium';
xlabel(ax4, 'Demand reactions');
ylabel(ax4, 'Medium component');

subplot(3,2,5);
heatmap(excess_redcost)
ax5 = gca;
ax5.XData = reaction_labels;
ax5.YData = medium_labels;
ax5.Title = 'Reduced cost in excess medium';
xlabel(ax5, 'Demand reactions');
ylabel(ax5, 'Medium component');

subplot(3,2,6);
heatmap(depletion_redcost)
ax6 = gca;
ax6.XData = reaction_labels;
ax6.YData = medium_labels;
ax6.Title = 'Reduced cost in depleted medium';
xlabel(ax6, 'Demand reactions');
ylabel(ax6, 'Medium component');

base = strcat('acetylationmodel_nucleus_', string(epsilon2));
fig_str = strcat(base, '.fig');
png_str = strcat(base, '.png');

saveas(fig, fig_str);
saveas(fig, png_str);
end