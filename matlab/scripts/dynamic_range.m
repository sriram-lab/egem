function [epsilon2] = dynamic_range(eps1, eps2, eps3, eps4, eps5, eps6, flag)

%% dynamic_range
    % Function returns the epsilon values that will be used for FBA
    % and FVA in metabolic_sensitivity.

    % INPUT:
        % eps1 to eps6: arrays containing values corresponding to fixed epsilon values
        % eps2: an array of epsilon values from 1E-6 to 1

    % OUTPUT:
        % epsilon2: an array that contains the epsilon values for
        % simultaneous optimization of metabolic model

if flag == 'excess'
    for i=1:6
        eval(["mx(i,:) = max(eps" i ".excess_flux);"]); 
        eval(["mn(i,:) = min(eps" i ".excess_flux);"]);
        dif(i,:) = mx(i,:) - mn(i,:);
    end
    
    % Get the epsilon2 that produced the highest dynamic range
    [~, idx] = max(dif);
    for i=1:length(idx)
        epsilon2(i,1) = eps2(idx(1,i));
    end
    
elseif flag == 'depletion'
    for i=1:6
        eval(["mx(i,:) = max(eps" i ".depletion_flux);"]);
        eval(["mn(i,:) = min(eps" i ".depletion_flux);"]);
        dif(i,:) = mx(i,:) - mn(i,:);
    end
    
    % Get the epsilon2 that produced the highest dynamic range
    [~, idx] = max(dif);
    for i=1:length(idx)
        epsilon2(i,1) = eps2(idx(1,i));
    end
end
end