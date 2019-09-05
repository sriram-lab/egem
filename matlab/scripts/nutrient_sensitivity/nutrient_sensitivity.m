% nutrient_sensitivity evaluates the impact of different nutrient
% conditions on the production of various histone markers
% @author: Scott Campit
function STRUCT = nutrient_sensitivity(model, medium, epsilon2, exp)

    load ./../../vars/nutrient_sensitivity
    reactions_of_interest = {'DM_KAC'; 'DM_KMe1'; 'DM_KMe2'; 'DM_KMe3'};
    GLUCOSE_UPTAKE_POS = find(ismember(model.rxns, 'EX_glc_e')); 
    tmp = model;
    
    for kappatype = 1:2
        constrained_model = medium_LB_constraints(tmp, medium);

        for component = 1:length(mediareactions1) % 50 medium components
            if kappatype == 1 
                excess_model = constrained_model;
                [~, pos]  = ismember(mediareactions1(:, 1), excess_model.rxns);
                
                if ismember(component, [1,4]) % glucose or glutamine
                    kappa = 3;
                    excess_model.lb(pos) = -media_exchange1(component)*kappa;
                else
                    kappa = 10;
                    excess_model.lb(pos) = -media_exchange1(component)*kappa;
                end

                switch exp
                    case {'sra', 'no_competition'}
                        [grate, flux] = singleReactionAnalysis(excess_model, ...
                            reactions_of_interest, epsilon2(:, 1));
                        
                        excess_grate(component, :) = grate;
                        excess_flux(component, :) = flux;
                        disp([component])
                                                
                    case 'competition'
                        [grate, flux] = competitiveFBA(excess_model, ...
                            reactions_of_interest, epsilon2(:, 1));
                        
                        excess_grate(component, :) = grate;
                        excess_flux(component, :) = flux;
                        disp(component)

                    case 'fva'
                        reaction_positions = find(ismember(excess_model.rxns, ...
                            reactions_of_interest));
                        excess_model.c(reaction_positions) = epsilon2(:, 1);
                        [~, maxFlux] = fluxVariability(excess_model, 100, ...
                            'max', reactions_of_interest);

                        excess_maxflux(component, :) = maxFlux;
                        disp(component)
                end

            elseif kappatype == 2 
                depletion_model = constrained_model;

                [~, pos] = ismember(mediareactions1(:, 1), depletion_model.rxns);
                if ismember(component, [2, 3, 5:19]) % trace nutrients
                    kappa = 0.0001;
                    depletion_model.lb(pos) = -media_exchange1(component)*kappa;
                else
                    kappa = 0.01;
                    depletion_model.lb(pos) = -media_exchange1(component)*kappa;
                end

                switch exp
                    case {'sra', 'no_competition'}
                        [grate, flux] = singleReactionAnalysis(depletion_model, ...
                            reactions_of_interest, epsilon2(:, 2));
                        
                        depletion_grate(component, :) = grate;
                        depletion_flux(component, :) = flux;
                        disp(component)

                    case 'competition'

                        [grate, flux] = competitiveFBA(depletion_model, ...
                            reactions_of_interest, epsilon2(:, 2));
                        
                        depletion_grate(component, :) = grate;
                        depletion_flux(component, :) = flux;
                        disp(component)

                    case 'fva'
                        reaction_positions = find(ismember(depletion_model.rxns, ...
                            reactions_of_interest));
                        depletion_model.c(reaction_positions) = epsilon2(:, 2);
                        [~, maxFlux] = fluxVariability(depletion_model, 100, ...
                            'max', reactions_of_interest);

                        depletion_maxflux(component, :) = maxFlux;
                        disp(component)
                end
            end
        end
    end

    % Initialize some structure variables
    STRUCT = struct('Name', exp);

    switch exp
        case {'sra', 'competition', 'no_competition'}
            medium_labels = mediareactions1(:, 2);
            reaction_labels = reactions_of_interest;
            fields = {...
                'excess_grate'; 'depletion_grate';...
                'excess_flux'; 'depletion_flux';...
                'medium_components'; 'reactions' ...
                };

            values = {...
                excess_grate; depletion_grate;...
                excess_flux; depletion_flux;...
                medium_labels; reaction_labels; ...
                };

            for i=1:length(fields)
                STRUCT.(fields{i}) = values{i};
            end

        case 'fva'
            medium_labels = mediareactions1(:,2);
            reaction_labels = reactions_of_interest;
            fields = {...
                'excess_maxflux'; 'depletion_maxflux'; ...
                'medium_components'; 'reactions'...
                };
            values = {...
                excess_maxflux; depletion_maxflux; ...
                medium_labels; reaction_labels...
                };
            for i=1:length(fields)
                STRUCT.(fields{i}) = values{i};
            end
    end
end
