function [rho, pvalue, R2, RMSE, MAE] ...
        = crossValidation(X, Y, kfold, trainSize)

for k = 1:kfold
        [Xtrain, Ytrain, Xtest, Ytest] = trainTestSplit(X, Y, trainSize);
        mdl = fitlm(Xtrain, Ytrain);
        YpredArr(k, :) = predict(mdl, Xtest);
        YtestArr(k, :) = Ytest;
end

finalYtest = reshape(YtestArr, [], 1);
finalYpred = reshape(YpredArr, [], 1);
[rho, pvalue] = corr(finalYtest, finalYpred);
RMSE = sqrt(mean(finalYtest - finalYpred).^2);
MAE = mean(abs(finalYtest - finalYpred));
R2 = 1 - (sum(finalYtest - finalYpred).^2 / ...
        sum(finalYtest - mean(finalYpred)).^2);

end