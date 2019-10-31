function [beta, costVector] = gradientDescent(X, Y, beta, alpha, iterations)
    costVector = zeros(iterations, size(Y, 2));
    for iter = 1:iterations
        prediction = ((X*beta - Y)' * X)';
        beta = beta - ((alpha/size(X, 1))*prediction);
        costVector(iter, :) = MeanSquareErr(X, Y, beta);
    end
end