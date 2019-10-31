function [Xtrain, Ytrain, Xtest, Ytest] = trainTestSplit(X, Y, trainSize)
    allVals = 1:size(Y, 1);
    allVals = transpose(allVals);
    vectorLength = round(trainSize*size(Y, 1));
    trainIndex = datasample(1:size(Y,1), vectorLength, 'Replace', false)';
    idx = ismember(allVals, trainIndex);
    testIndex = allVals(~idx);
        
    Xtest = X(testIndex, :);
    Xtrain = X(trainIndex, :);
    Ytest = Y(testIndex, 1);
    Ytrain = Y(trainIndex, 1);
end