%% @author: Scott Campit
function STRUCT = metabolic_sensitivity(model, reactions_of_interest,...
    compartment, epsilon2, scaling, exp, medium, fva_grate)
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

if (~exist('reactions_of_interest', 'var')) || (isempty('reactions_of_interest'))
    load('./../vars/metabolites.mat')
end


% Load substrate uptake rates, medium components, reactions of interest
load ./../vars/supplementary_software_code media_exchange1
load('./../vars/mediareactions1.mat') % Medium components
metabolites = reactions_of_interest;

% Initialize parameters needed
minfluxflag = 0;
GLUCOSE_UPTAKE_POS = find(ismember(model.rxns, 'EX_glc(e)'));  % glucose uptake reaction in eGEM model
BIOMASS_OBJ_POS = find(ismember(model.rxns, 'biomass_objective')); % biomass rxn position in eGEM model
rxnName = metabolites(:, 1); % reaction positions of interest
rxnName = char(rxnName);

% SRA: Same epsilon value for all reactions. 
if isequal(size(epsilon2),[1,1])
    epsilon2 = repmat(epsilon2, 20, 2);
end

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
                
                % Single reaction optimization and capturing results for
                % the optimal FBA conditions
                case {'sra', 'no_competition'}
                    
                    % Add demand reactions from the metabolite list to the metabolic model
                    for rxn = 1:length(metabolites(:,1)) 

                        % Create the demand reactions dynamically
                        tmp_met = char(metabolites(rxn, 2)); % example: accoa
                        met = [tmp_met '[' compartment ']'];
                        tmprxn = [tmp_met '[' compartment '] -> ']; % example: accoa[n] ->
                        tmpname = char(metabolites(rxn, 1)); % example: DM_accoa
                        
                        % Add the new reaction to the excess and depletion models
                        excess_model = addReaction(excess_model, tmpname,...
                            'reactionFormula', tmprxn);

                        % Calculate growth rates, reduced costs, and shadow prices
                        [soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(excess_model,...
                            BIOMASS_OBJ_POS, met, [], [], [], 1, exp);
                        excess_grate(component, rxn) = grate;
                        excess_grate_sp(component, rxn) = grate_sp;
                        excess_grate_rc(component, rxn) = grate_rc;

                        %% Obtain flux values when using epsilon2 as the objective coefficient for the reaction of interest.
                        model3 = excess_model;

                        % Set objective of the reaction of interest to epsilon2
                        rxnpos = [find(ismember(model3.rxns, tmpname))];
                        epsilon2_excess = epsilon2(:, 1);
                        model3.c(rxnpos) = epsilon2_excess(rxn, 1);

                        % Solve for metabolic fluxes
                        [soln, flux, flux_sp, flux_rc, ~, ~] =...
                            calc_metabolic_metrics(model3,...
                            rxnpos, met, [], [], [], epsilon2_excess(rxn, 1), exp);
                        excess_flux(component, rxn) = flux;
                        excess_flux_sp(component,rxn) = flux_sp;
                        excess_flux_rc(component,rxn) = flux_rc;

                        % Reset models
                        excess_model = removeRxns(excess_model, tmprxn);
                        disp([component, rxn])
                    end
     
                % FBA optimization for all reactions simultaneously
                case 'competition'
                    
                    % Get all the reactions you are interested in
                    mets = cellstr(metabolites(:, 2)); 
                    for i=1:length(mets)
                        met{i} = [char(mets(i)) '[' compartment ']'];
                    end
                    met = string(met);

                    % Get growth rate and other metrics from the function
                    [excess_soln, grate, grate_sp, grate_rc, ~, ~]...
                        = calc_metabolic_metrics(excess_model,...
                            BIOMASS_OBJ_POS, met, [], [], [], 1, exp);

                    % Stuff you want
                    excess_grate(component,:) = grate;
                    excess_grate_sp(component,:) = grate_sp;
                    excess_grate_rc(component,:) = grate_rc;

                    %% Simultaneously solve several reactions and get the flux
                    model3 = excess_model;
                    rxnpos = [find(ismember(model3.rxns, rxnName))];

                    % Get metabolic flux and other metrics from the function
                    [excess_soln, flux, flux_sp, flux_rc, ~, ~] =...
                        calc_metabolic_metrics(model3, rxnpos, met, [],...
                        [], [], epsilon2(:, 1), exp);

                    % Stuff you want
                    excess_flux(component,:) = flux;
                    excess_flux_sp(component,:) = flux_sp;
                    excess_flux_rc(component,:) = flux_rc;

                    disp(component)

                % FVA optimization for all reactions simultaneously  
                case 'fva'

                    % Not interested in growth rates - just get the fluxes
                    model3 = excess_model;
                    rxnpos = [find(ismember(model3.rxns, rxnName))];
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

                % Single reaction optimization and capturing results for
                % the optimal FBA conditions
                case {'sra', 'no_competition'}
                    
                    % Add demand reactions from the metabolite list to the metabolic model
                    for rxn = 1:length(metabolites(:,1)) 

                        % Create the demand reactions dynamically
                        tmp_met = char(metabolites(rxn, 2)); % example: accoa
                        met = [tmp_met '[' compartment ']'];
                        tmprxn = [tmp_met '[' compartment '] -> ']; % example: accoa[n] ->
                        tmpname = char(metabolites(rxn, 1)); % example: DM_accoa
                        
                        % Add the new reaction to the excess and depletion models
                        depletion_model = addReaction(depletion_model, tmpname,...
                            'reactionFormula', tmprxn);

                        % Calculate growth rates, reduced costs, and shadow prices
                        [soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(depletion_model,...
                            BIOMASS_OBJ_POS, met, [], [], [], 1, exp);
                        depletion_grate(component, rxn) = grate;
                        depletion_grate_sp(component, rxn) = grate_sp;
                        depletion_grate_rc(component, rxn) = grate_rc;

                        %% Obtain flux values when using epsilon2 as the objective coefficient for the reaction of interest.
                        model3 = depletion_model;

                        % Set objective of the reaction of interest to epsilon2
                        rxnpos = [find(ismember(model3.rxns, tmpname))];
                        epsilon2_depletion = epsilon2(:, 2);
                        model3.c(rxnpos) = epsilon2_depletion(rxn);

                        % Solve for metabolic fluxes
                        [soln, flux, flux_sp, flux_rc, ~, ~] = calc_metabolic_metrics(model3,...
                            rxnpos, met, [], [], [], epsilon2_depletion(rxn), exp);
                        depletion_flux(component, rxn) = flux;
                        depletion_flux_sp(component,rxn) = flux_sp;
                        depletion_flux_rc(component,rxn) = flux_rc;

                        % Reset models
                        depletion_model = removeRxns(depletion_model, tmprxn);
                        disp([component, rxn])
                    end

                % FBA optimization for all reactions simultaneously
                case 'competition'

                    mets = cellstr(metabolites(:, 2)); 
                    for i=1:length(mets)
                        met{i} = [char(mets(i)) '[' compartment ']'];
                    end
                    met = string(met);

                    [depletion_soln, grate, grate_sp, grate_rc, ~, ~] = calc_metabolic_metrics(depletion_model,...
                            BIOMASS_OBJ_POS, met, [], [], [], 1, exp);

                    depletion_grate(component, :) = grate;
                    depletion_grate_sp(component, :) = grate_sp;
                    depletion_grate_rc(component, :) = grate_rc;

                    %% Simultaneously solve several reactions and get the flux
                    model3 = depletion_model;
                    rxnpos = [find(ismember(model3.rxns, rxnName))];

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
                    rxnpos = [find(ismember(model3.rxns, rxnName))];
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

% Initialize some structure variables
STRUCT = struct('Name', exp);

%% Saving data structures, generating figures, etc.
switch exp
    case {'sra', 'competition', 'no_competition'}
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3)';
        fields = {...
            'excess_grate'; 'depletion_grate';...
            'excess_grate_sp'; 'depletion_grate_sp';...
            'excess_grate_rc'; 'depletion_grate_rc';...
            'excess_flux'; 'depletion_flux';...
            'excess_flux_sp'; 'depletion_flux_sp';...
            'excess_flux_rc'; 'depletion_flux_rc';...
            'medium_components'; 'reactions' ...
            };
        
        values = {...
            excess_grate; depletion_grate;...
            excess_grate_sp; depletion_grate_sp;...
            excess_grate_rc; depletion_grate_rc;...
            excess_flux; depletion_flux;...
            excess_flux_sp; depletion_flux_sp;...
            excess_flux_rc; depletion_flux_rc;...
            medium_labels; reaction_labels; ...
            };
        
        for i=1:length(fields)
            STRUCT.(fields{i}) = values{i};
        end

        base = strcat('./../tables/eGEM_', exp);   
        table_str = strcat(base, '.xlsx');
        rowname = medium_labels;
        colname = reaction_labels;
        %save_xl18(table_str, rowname, colname, sra, epsilon2, exp)

    case 'fva'
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);
        fields = {...
            'excess_maxflux'; 'excess_minflux';...
            'depletion_maxflux'; 'depletion_minflux';...
            'medium_components'; 'reactions'...
            };
        values = {...
            excess_maxflux; excess_minflux;...
            depletion_maxflux; depletion_minflux;...
            medium_labels; reaction_labels...
            };
        for i=1:length(fields)
            STRUCT.(fields{i}) = values{i};
        end
end

end
