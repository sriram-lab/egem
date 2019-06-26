function [epsilon2] = dynamic_range(arr1, arr2, arr3, arr4, arr5, arr6, eps2)
% dynamic_range returns the epsilon values that will be used for FBA
% and FVA in metabolic_sensitivity.

% INPUT:
    % arr1 to arr6: arrays containing values corresponding to fixed epsilon values
    % eps2: an array of epsilon values from 1E-6 to 1
    
% OUTPUT:
    % epsilon2: an array that contains the epsilon values for
    % simultaneous optimization of metabolic model

% For each array, get the dynamic range, stored in `dif` variable
for i=1:6
    eval(["mx(i,:) = max(arr" i ");"]); %(1x20) array
    eval(["mn(i,:) = min(arr" i ");"]);
    dif(i,:) = mx(i,:) - mn(i,:);
end

% Get the epsilon2 that produced the highest dynamic range
[~, idx] = max(dif);
for i=1:length(idx)
    epsilon2(i,1) = eps2(idx(1,i));
end
end