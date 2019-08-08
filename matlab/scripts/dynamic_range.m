% `dynamic_range`: identifies the reactions that have the largest dynamic
% flux range by adding or removing medium components to the medium
% condition
% @author: Scott Campit
function epsilon2 = dynamic_range(eps1, eps2, eps3, eps4, eps5, eps6, eps7,...
    type)
    
    function [dynamic_range, maximum, minimum] = compute_dynamic_range(eps1, ...
            eps2, eps3, eps4, eps5, eps6, eps7, medium_type)
        for structure_number = 1:7
            eval(strcat("maximum(structure_number, :) = max(eps",...
                string(structure_number), ".", medium_type, "_flux);")); 
            eval(strcat("minimum(structure_number, :) = min(eps", ...
                string(structure_number), ".", medium_type, "_flux);")); 
            dynamic_range(structure_number, :) = ...
                maximum(structure_number, :) - minimum(structure_number, :);
        end
    end
    
    function objCoefs = getObjCoefs(values, objective_coefficients)
        [~, index] = max(values);
        for structure_number = 1:length(index)
            objCoefs(structure_number, 1) = ...
                objective_coefficients(index(1, structure_number));
        end
    end
    
    objective_coefficients = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1];
    switch type
        case 'dynamic'
            [excess_dynamic_range, ~, ~] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'excess');
            excess_objCoefs = getObjCoefs(excess_dynamic_range, ...
                objective_coefficients);
            
            [depletion_dynamic_range, ~, ~] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'depletion');
            depletion_objCoefs = getObjCoefs(depletion_dynamic_range, ...
                objective_coefficients);
            
            epsilon2 = [excess_objCoefs, depletion_objCoefs];

        case 'max'
            [~, excess_maximum, ~] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'excess');
            excess_objCoefs = getObjCoefs(excess_maximum, ...
                objective_coefficients);
            
            [~, depletion_maximum, ~] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'excess');
            depletion_objCoefs = getObjCoefs(depletion_maximum, ...
                objective_coefficients);
            
            epsilon2 = [excess_objCoefs, depletion_objCoefs];

        case 'min'
            [~, ~, excess_minimum] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'excess');
            excess_objCoefs = getObjCoefs(excess_minimum, ...
                objective_coefficients);
            
            [~, ~, depletion_minimum] = compute_dynamic_range(eps1, eps2, ...
                eps3, eps4, eps5, eps6, eps7, 'excess');
            depletion_objCoefs = getObjCoefs(depletion_minimum, ...
                objective_coefficients);
            
            epsilon2 = [excess_objCoefs, depletion_objCoefs];
    end
end