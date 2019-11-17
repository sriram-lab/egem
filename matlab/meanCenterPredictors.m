function [centered] = meanCenterPredictors(array)
    for col = 1:size(array, 2)
       centered(:, col) = normalize(array(:, col));
    end
end
