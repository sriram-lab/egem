function [betaHat, Yhat, error, CovResidual] = multiVarRegress(X, Y)
    betaHat = inv(X'*X)*X'*Y;
    Yhat = X*betaHat;
    error = Y - Yhat;
    
    b = betaHat(:);
    
    %for i = 1:length(b)
    %    for j = 1:length(b)
    %        CovB(i,j) = ((b(i)-mean(b))*(b(j)-mean(b)'));
    %    end
    %end
    CovResidual = cov(error,1);
end