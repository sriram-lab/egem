function [CVResidual] = crossValidation(X, Y, ntimes, maxSplit)
    for i = 1:ntimes
        index = randomperm(round(size(X, 1)*maxSplit));
        
        Xtest = X(~index, :);
        Xtrain = X(index, :);
        Ytest = Y(~index 1);
        Ytrain = Y(index, 1);
        
        [beta, betaInt, r, rint, stats] = regress(Ytrain, Xtrain);
        YCV = Xtest*beta;
        CVResidual = YCV - Ytest;
    end
end