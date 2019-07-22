%% dynamic_range
function epsilon2 = dynamic_range(eps1, eps2, eps3, eps4, eps5, eps6, eps7, type)
    % Dynamic Range
    % Function returns the epsilon values that will be used for FBA
    % and FVA in metabolic_sensitivity.

    % INPUT:
        % eps1 to eps6: arrays containing values corresponding to fixed epsilon values
        % eps2: an array of epsilon values from 1E-6 to 1

    % OUTPUT:
        % epsilon2: an array that contains the epsilon values for
        % simultaneous optimization of metabolic model
switch type
    case 'dynamic'
        % Dynamic range for metabolic flux in excess medium conditions
        for i=1:7
            eval(strcat("mx(i,:) = max(eps", string(i), ".excess_flux);")); 
            eval(strcat("mn(i,:) = min(eps", string(i), ".excess_flux);"));
            dif(i,:) = mx(i,:) - mn(i,:);
        end

        [~, idx] = max(dif);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            excess_epsilon2(i,1) = fill_val(idx(1,i));
        end

        % Dynamic range for metabolic flux in depleted medium conditions
        for i=1:7
            eval(strcat("mx(i,:) = max(eps", string(i), ".depletion_flux);")); 
            eval(strcat("mn(i,:) = min(eps", string(i), ".depletion_flux);"));
            dif(i,:) = mx(i,:) - mn(i,:);
        end

        [~, idx] = max(dif);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            depletion_epsilon2(i,1) = fill_val(idx(1,i));
        end

        % Concatenate
        epsilon2 = [excess_epsilon2, depletion_epsilon2];
        
    case 'max'
        % Dynamic range for metabolic flux in excess medium conditions
        for i=1:7
            eval(strcat("mx(i,:) = max(eps", string(i), ".excess_flux);")); 
        end

        [~, idx] = max(mx);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            excess_epsilon2(i,1) = fill_val(idx(1,i));
        end
        
        % Dynamic range for metabolic flux in depleted medium conditions
        for i=1:7
            eval(strcat("mx(i,:) = max(eps", string(i), ".depletion_flux);")); 
        end

        [~, idx] = max(mx);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            depletion_epsilon2(i,1) = fill_val(idx(1,i));
        end
        
        % Concatenate
        epsilon2 = [excess_epsilon2, depletion_epsilon2];
        
    case 'min'
        % Dynamic range for metabolic flux in excess medium conditions
        for i=1:7
            eval(strcat("mn(i,:) = min(eps", string(i), ".excess_flux);"));
        end

        [~, idx] = max(mn);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            excess_epsilon2(i,1) = fill_val(idx(1,i));
        end
        
        % Dynamic range for metabolic flux in depleted medium conditions
        for i=1:7
            eval(strcat("mn(i,:) = min(eps", string(i), ".depletion_flux);"));
        end

        [~, idx] = max(mn);
        fill_val = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
        for i=1:length(idx)
            depletion_epsilon2(i,1) = fill_val(idx(1,i));
        end
        
        % Concatenate
        epsilon2 = [excess_epsilon2, depletion_epsilon2];
end
end