%% @author: Scott Campit
function [sra, fba, fva, excess_soln, depletion_soln] = metabolic_sensitivity(model, compartment,...
    epsilon2, scaling, exp, medium, fva_grate)
% metabolic_sensitivity.m displays the values corresponding to several demand
% reactions and excess/depletion of a specific medium component.

% INPUTS:
    % model: genome scale metabolic model
    % rxnpos: list of reactions that will be used in the analysis
    % epsilon2: a weight that will be used for the secondary obcomponentective coefficient
    % scaling: if the data should be scaled/normalized
    % exp: this flag can take 3 possible values:
        % 'sra': optimizes single reactions with a fixed epsilon value
        % 'fba': optimizes multiple reactions using fba
        % 'fva': sets the biomass to be at 100% and gets the maximum
        % metabolic fluxes associated

% OUTPUTS:
    % excess_flux: metabolic flux corresponding to excess medium
    % depletion_flux: ''' depletion medium
    % excess_redcost: reduced cost corresponding to excess medium
    % depletion_redcost: ''' depletion medium
    % excess_shadow: shadow costs corresponding to excess medium
    % depletion_shadow: ''' depletion medium

%% metabolic_sensitivity.m

% Default is RPMI medium
%if (~exist('medium','var')) || (isempty(medium))
%    medium = 'RPMI';
%end

% Load substrate uptake rates, medium components, reactions of interest
load ./../vars/supplementary_software_code media_exchange1
var = {'./../vars/metabolites.mat', './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk});
end

% Initialize parameters needed
minfluxflag = 0;
posgluc = find(ismember(model.rxns, 'EX_glc(e)'));  % glucose uptake reaction in eGEM model
biomassobjpos = find(ismember(model.rxns, 'biomass_objective')); % biomass rxn position in eGEM model
rxnname = metabolites(:, 1); % reaction positions of interest
rxnname(5) = {'EX_KAC'}; % I need to rename this reaction in the metabolites folder
rxnname = char(rxnname);

% Set the scaling proportion for RPMI substrate uptake rates to
% a fixed amount. weight = 10 is excess, and weight = 0.01 is
% depletion
for kappatype = 1:2
    
    % Create a temporary metabolic model that will have substrate uptake
    % rates fitted on
    tmp = model;
    
    % Make the methylation exchange reaction have a fixed LB of
    % -0.5 to be non-limiting
    [~, pos] = ismember({'EX_met_L(e)'}, tmp.rxns);
    tmp.lb(pos) = -0.5;    
    
    % Set the substrate uptake rates based on medium
    % formulation (taken from medium sources). 
    tmp = media(tmp, medium);
    
    % For each medium component, set the substrate uptake rate. 
    for component = 1:length(mediareactions1(:,1)) % 50 medium components
        %% Excess medium conditions
        if kappatype == 1 % glucose or glutamine
            excess_model = tmp;
            weight  = 10;
            
            [~, pos]  = ismember(mediareactions1(component,1), excess_model.rxns);
            if ismember(component, [1;4])
                kappa = 3;
                excess_model.lb(pos) = -media_exchange1(component,1)*kappa;
            else
                kappa = weight;
                excess_model.lb(pos) = -media_exchange1(component,1)*kappa;
            end

            switch exp
            
                % Single reaction optimization
                case 'sra'
                    fba = 0;
                    fva = 0;

                    % Add demand reactions from the metabolite list to the metabolic model
                    for rxn = 1:length(metabolites(:,1)) % 20 reactions of interest

                        % Create the demand reactions dynamically
                        tmp_met = char(metabolites(rxn, 2)); % example: accoa
                        tmprxn = [tmp_met '[' compartment '] -> ']; % example: accoa[n] ->
                        tmpname = char(metabolites(rxn, 1)); % example: DM_accoa

                        % Add the new reaction to the excess and depletion models
                        excess_model = addReaction(excess_model, tmpname,...
                            'reactionFormula', tmprxn);

                        % Calculate growth rates, reduced costs, and shadow prices
                        [soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(excess_model,...
                            biomassobjpos, tmp_met, [], [], [], 1, exp);
                        excess_grate(component, rxn) = grate;
                        excess_grate_sp(component, rxn) = grate_sp;
                        excess_grate_rc(component, rxn) = grate_rc;

                        %% Obtain flux values when using epsilon2 as the objective coefficient for the reaction of interest.
                        model3 = excess_model;

                        % Set objective of the reaction of interest to epsilon2
                        rxnpos = [find(ismember(model3.rxns, tmpname))];
                        model3.c(rxnpos) = epsilon2;

                        % Solve for metabolic fluxes
                        [flux, flux_sp, flux_rc, ~, ~] = calc_metabolic_metrics(model3,...
                            rxnpos, tmp_met, [], [], [], epsilon2, exp);
                        excess_flux(component, rxn) = flux;
                        excess_flux_sp(component,rxn) = flux_sp;
                        excess_flux_rc(component,rxn) = flux_rc;

                        % Reset models
                        excess_model = removeRxns(excess_model, tmprxn);
                        disp([component, rxn])
                    end
                    
                % FBA optimization for all reactions simultaneously
                case 'fba'
                    sra=0;
                    fva=0;

                    % Get all the reactions you are interested in
                    mets = cellstr(metabolites(:, 2)); 
                    for i=1:length(mets)
                        met{i} = [char(mets(i)) '[' compartment ']'];
                    end
                    met = string(met);

                    % Get growth rate and other metrics from the function
                    [grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(excess_model,...
                            biomassobjpos, met, [], [], [], 1, exp);

                    % Stuff you want
                    excess_grate(component,:) = grate';
                    excess_grate_sp(component,:) = grate_sp';
                    excess_grate_rc(component,:) = grate_rc';

                    %% Simultaneously solve several reactions and get the flux
                    model3 = excess_model;
                    rxnpos = [find(ismember(model3.rxns, rxnname))];

                    % Get metabolic flux and other metrics from the function
                    [excess_soln, flux, flux_sp, flux_rc, ~, ~] = calc_metabolic_metrics(model3,...
                            rxnpos, met, [], [], [], epsilon2(:, 1), exp);

                    % Stuff you want
                    excess_flux(component,:) = flux';
                    excess_flux_sp(component,:) = flux_sp';
                    excess_flux_rc(component,:) = flux_rc';

                    disp(component)

                % FVA optimization for all reactions simultaneously  
                case 'fva'
                    sra=0;
                    fba=0;
                    excess_soln = 0;
                    depletion_soln = 0;

                    % Not interested in growth rates - just get the fluxes
                    model3 = excess_model;
                    rxnpos = [find(ismember(model3.rxns, rxnname))];
                    rxn_nam = model3.rxns(rxnpos);

                    % epsilon values for excess
                    model3.c(rxnpos) = epsilon2(:, 1);

                    % Run FVA for methylation and acetylation reactions.
                    [~, ~, ~, ~, maxflux, minflux] = calc_metabolic_metrics(model3,...
                            rxnpos, [], fva_grate, 'max', rxn_nam, epsilon2(:, 1), exp);

                    excess_maxflux(component, :) = maxflux;
                    excess_minflux(component, :) = minflux;
                    disp(component)
            end

            %% Depleted medium conditions
            elseif kappatype == 2 % trace elements
                depletion_model = tmp;
                weight = 0.01;

                [~, pos]  = ismember(mediareactions1(component,1), depletion_model.rxns);

                if ismember(component,[2,3,5:19])
                    kappa = weight/100;
                    depletion_model.lb(pos) = -media_exchange1(component,1)*kappa;
                else
                    kappa = weight;
                    depletion_model.lb(pos) = -media_exchange1(component,1)*kappa;
                end

                switch exp

                    % Single reaction optimization
                    case 'sra'
                        % Add demand reactions from the metabolite list to the metabolic model
                        for rxn = 1:length(metabolites(:,1)) % 20 reactions of interest

                            % Create the demand reactions dynamically
                            tmp_met = char(metabolites(rxn, 2)); % example: accoa
                            tmprxn = [tmp_met '[' compartment '] -> ']; % example: accoa[n] ->
                            tmpname = char(metabolites(rxn, 1)); % example: DM_accoa

                            % Add the new reaction to the excess and depletion models
                            depletion_model = addReaction(depletion_model, tmpname,...
                                'reactionFormula', tmprxn);

                            % Calculate growth rates, reduced costs, and shadow prices
                            [soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(depletion_model,...
                                biomassobjpos, tmp_met, [], [], [], 1, exp);
                            depletion_grate(component, rxn) = grate;
                            depletion_grate_sp(component,rxn) = grate_sp;
                            depletion_grate_rc(component,rxn) = grate_rc;

                            %% Obtain flux values when using epsilon2 as the objective coefficient for the reaction of interest.
                            model3 = depletion_model;

                            % Set objective of the reaction of interest to epsilon2
                            rxnpos = [find(ismember(model3.rxns, tmpname))];
                            model3.c(rxnpos) = epsilon2;

                            % Solve for metabolic fluxes
                            [soln, flux, flux_sp, flux_rc, ~, ~] = calc_metabolic_metrics(model3,...
                                rxnpos, tmp_met, [], [], [], epsilon2, exp);
                            depletion_flux(component, rxn) = flux;
                            depletion_flux_sp(component,rxn) = flux_sp;
                            depletion_flux_rc(component,rxn) = flux_rc;

                            % Reset models
                            depletion_model = removeRxns(depletion_model, tmprxn);
                            disp([component, rxn])
                        end

                    % FBA optimization for all reactions simultaneously
                    case 'fba'

                        mets = cellstr(metabolites(:, 2)); 
                        for i=1:length(mets)
                            met{i} = [char(mets(i)) '[' compartment ']'];
                        end
                        met = string(met);

                        [depletion_soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(depletion_model,...
                                biomassobjpos, met, [], [], [], 1, exp);

                        depletion_grate(component, :) = grate;
                        depletion_grate_sp(component, :) = grate_sp;
                        depletion_grate_rc(component, :) = grate_rc;

                        %% Simultaneously solve several reactions and get the flux
                        model3 = depletion_model;
                        rxnpos = [find(ismember(model3.rxns, rxnname))];

                        [soln, flux, flux_sp, flux_rc, ~, ~] = calc_metabolic_metrics(model3,...
                                rxnpos, met, [], [], [], epsilon2(:, 2), exp);

                        depletion_flux(component, :) = flux;
                        depletion_flux_sp(component, :) = flux_sp;
                        depletion_flux_rc(component, :) = flux_rc;

                        disp(component)

                    % FVA optimization for all reactions simultaneously  
                    case 'fva'
                        % Not interested in growth rates - just get the fluxes
                        model3 = depletion_model;
                        rxnpos = [find(ismember(model3.rxns, rxnname))];
                        rxn_nam = model3.rxns(rxnpos);

                        % epsilon values for depletion
                        model3.c(rxnpos) = epsilon2(:, 2);

                        % Run FVA for methylation and acetylation reactions.
                        [~, ~, ~, ~, maxflux, minflux] = calc_metabolic_metrics(model3,...
                                rxnpos, [], fva_grate, 'max', rxn_nam, epsilon2(:, 2), exp);
                        depletion_maxflux(component, :) = maxflux;
                        depletion_minflux(component, :) = minflux;
                        disp(component)
                end
            end
        end
end

%% Package results into prettified structures
if exp == 'sra'
    medium_labels = mediareactions1(:,2);
    reaction_labels = metabolites(:,3);
    fields = {...
        'excess_grate'; 'depletion_grate';...
        'excess_grate_sp'; 'depletion_grate_sp';...
        'excess_grate_rc'; 'depletion_grate_rc';...
        'excess_flux'; 'depletion_flux';...
        'excess_flux_sp'; 'depletion_flux_sp';...
        'excess_flux_rc'; 'depletion_flux_rc';...
        %'medium_components'; 'reactions'
        };
    values = {...
        excess_grate; depletion_grate;...
        excess_grate_sp; depletion_grate_sp;...
        excess_grate_rc; depletion_grate_rc;...
        excess_flux; depletion_flux;...
        excess_flux_sp; depletion_flux_sp;...
        excess_flux_rc; depletion_flux_rc;...
        %medium_labels; reaction_labels...
        };
    for i=1:length(fields)
        fba.(fields{i}) = values{i};
    end
    
elseif exp == 'fba'
    medium_labels = mediareactions1(:,2);
    reaction_labels = metabolites(:,3);
    fields = {...
        'excess_grate'; 'depletion_grate';...
        'excess_grate_sp'; 'depletion_grate_sp';...
        'excess_grate_rc'; 'depletion_grate_rc';...
        'excess_flux'; 'depletion_flux';...
        'excess_flux_sp'; 'depletion_flux_sp';...
        'excess_flux_rc'; 'depletion_flux_rc';...
        %'medium_components'; 'reactions'...
        };
    values = {...
        excess_grate; depletion_grate;...
        excess_grate_sp; depletion_grate_sp;...
        excess_grate_rc; depletion_grate_rc;...
        excess_flux; depletion_flux;...
        excess_flux_sp; depletion_flux_sp;...
        excess_flux_rc; depletion_flux_rc;...
        %medium_labels; reaction_labels...
        };
    for i=1:length(fields)
        fba.(fields{i}) = values{i};
    end
    
elseif exp == 'fva'
    medium_labels = mediareactions1(:,2);
    reaction_labels = metabolites(:,3);
    fields = {...
        'excess_max_flux'; 'excess_min_flux';...
        'depletion_max_flux'; 'depletion_min_flux';...
        %'medium_components'; 'reactions'...
        };
    values = {...
        excess_maxflux; excess_minflux;...
        depletion_maxflux; depletion_minflux;...
        %medium_labels; reaction_labels...
        };
    for i=1:length(fields)
        fva.(fields{i}) = values{i};
    end
end

% %% Scaling the data for visualization purposes
if scaling == 'zscore'
    if exp == 'sra'
        sra = structfun(@zscore, sra, 'UniformOutput', false);
    elseif exp == 'fba'
        fba = structfun(@zscore, fba, 'UniformOutput', false);
    elseif exp == 'fva'
        fva = structfun(@zscore, fva, 'UniformOutput', false);
    end
end
 
% %% Replace 0 with NaN
% %excess_flux(excess_flux==0) = NaN;
% %excess_redcost(excess_redcost==0) = NaN;
% %excess_shadow(excess_shadow==0) = NaN;
% 
% %depletion_flux(depletion_flux==0) = NaN;
% %depletion_redcost(depletion_redcost==0) = NaN;
% %depletion_shadow(depletion_shadow==0) = NaN;
 
%% Heatmap for FVA

switch exp
    % Case 1: Use the most dynamic range of metabolic fluxes
    case 'fva'  
        % Prepare figure labels and variables
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);

        fig = figure;
        subplot(1,2,1);
        heatmap(excess_maxflux)
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(depletion_maxflux)
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        % Save figure
        
        fig_str = strcat(string(medium), '.fig');
        saveas(fig, fig_str);
end
% 
% 
%% Heatmap figures
% Still needs work:
    % Make `excess` green and `depletion` red
    % Color gradient: metabolic flux > reduced costs > shadow prices
    % Incorporate plotly
    % Is there a way to dynamically size .pngs?

switch exp
    % Case 1: Use the most dynamic range of metabolic fluxes
    case 'fba'  
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
        heatmap(excess_flux_sp)
        ax3 = gca;
        ax3.XData = reaction_labels;
        ax3.YData = medium_labels;
        ax3.Title = 'Shadow price in excess medium';
        xlabel(ax3, 'Demand reactions');
        ylabel(ax3, 'Medium component');

        subplot(2,3,5);
        heatmap(depletion_flux_sp)
        ax4 = gca;
        ax4.XData = reaction_labels;
        ax4.YData = medium_labels;
        ax4.Title = 'Shadow price in depleted medium';
        xlabel(ax4, 'Demand reactions');
        ylabel(ax4, 'Medium component');

        subplot(2,3,3);
        heatmap(excess_flux_rc)
        ax5 = gca;
        ax5.XData = reaction_labels;
        ax5.YData = medium_labels;
        ax5.Title = 'Reduced cost in excess medium';
        xlabel(ax5, 'Demand reactions');
        ylabel(ax5, 'Medium component');

        subplot(2,3,6);
        heatmap(depletion_flux_rc)
        ax6 = gca;
        ax6.XData = reaction_labels;
        ax6.YData = medium_labels;
        ax6.Title = 'Reduced cost in depleted medium';
        xlabel(ax6, 'Demand reactions');
        ylabel(ax6, 'Medium component');
end
% 
%         base = strcat('./../figures/new-model/eGEMn_', string(epsilon2));
%         fig1_str = strcat(base, '.fig');
% 
%         saveas(fig1, fig1_str);
% 
%         % Create a heatmap for the growth rates in excess and depleted medium
%         fig2 = figure;
% 
%         dat = [excess_xgrate; depletion_xgrate];
%         dat = dat';
% 
%         heatmap(dat)
%         ax1 = gca;
%         %ax1.XData = ['Excess Medium'; 'Depleted Medium'];
%         ax1.YData = medium_labels;
%         ax1.Title = 'Growth rate in excess medium';
%         xlabel(ax1, 'Growth Rates');
%         ylabel(ax1, 'Medium component');
% 
%         subplot(1,2,1);
%         heatmap(excess_xgrate)
%         ax1 = gca;
%         ax1.XData = 'Growth Rate';
%         ax1.YData = medium_labels;
%         ax1.Title = 'Growth rate in excess medium';
%         xlabel(ax1, 'Demand reactions');
%         ylabel(ax1, 'Medium component');
% 
%         subplot(1,2,2);
%         heatmap(depletion_xgrate)
%         ax2 = gca;
%         ax2.XData = 'Growth Rate';
%         ax2.YData = medium_labels;
%         ax2.Title = 'Growth rate in depleted medium';
%         xlabel(ax2, 'Demand reactions');
%         ylabel(ax2, 'Medium component');
% 
%         base = strcat('./../figures/new-model/eGEMn_grate_', string(epsilon2));
%         fig2_str = strcat(base, '.fig');
%         saveas(fig2, fig2_str);
% end
% %clear
end
