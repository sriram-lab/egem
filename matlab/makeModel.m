function [model] = makeModel(data)
    % Make data arrays and compute relevant stuff
    array = table2array(data);
    X = meanCenterPredictors(array);
    Y = X(:, 2);
    X(:, 2) = [];
    trueBeta = X\Y;
    model = fitlm(X, Y);
end