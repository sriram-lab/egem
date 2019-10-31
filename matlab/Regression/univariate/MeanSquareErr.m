function [cost] = MeanSquareErr(X, Y, beta)
    prediction = X*beta;
    cost = 1/(2*size(X,1))*sum(prediction-Y).^2;
end