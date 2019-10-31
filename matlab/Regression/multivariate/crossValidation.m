function [CVResidual] = crossValidation(X, Y, ntimes, maxSplit)
    for i = 1:ntimes
        [Xtrain, Ytrain, Xtest, Ytest] = trainTestSplit(X, Y, 0.8);
        [beta, betaInt, r, rint, stats] = regress(Ytrain, Xtrain);
        YCV = Xtest*beta;
        CVResidual = YCV - Ytest;
    end
end