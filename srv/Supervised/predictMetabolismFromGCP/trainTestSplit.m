function data = trainTestSplit(X, Y, trainingSize, randomState)
%% TRAINTESTSPLIT trainTestSplit
% |trainTestSplit| separates the the predictors ( |X| ) and response ( |Y| ) 
% variables into a training and test set. 
% 
% *INPUTS:*
% 
% |X:|                                                  A numeric array of n 
% observations by p predictor variables.
% 
% |Y:|                                                  A numeric or boolean 
% array of n observations by 1 response variable.
% 
% |trainingSize (optional):|   A scalar value ranging from [0, 1] denoting the 
% size of the dataset for training. The default value is 80%.
% 
% |randomState (optional):|     A boolean / logical value that determines whether 
% should be shuffled using a fixed seed value. The default value is |true|.
% 
% *OUTPUTS:*
% 
% |Xtrain|:         An array that contains p predictors and n * trainingSize 
% observations for training the model.
% 
% |Xtest|:           An array that contains p predictors and n * (1 - trainingSize) 
% observations for testing the model.
% 
% |Ytrain|:         A vector containing n * trainingSize observations for training 
% the model.
% 
% |Ytest|:           A vector containing n * (1 - trainingSize) observations 
% for testing the model.
    
    if nargin < 3
        trainingSize = 0.8;
    end
    
    if nargin < 4
        % For reproducibility - if you want randomly shuffled data, turn this
        % off.
        rng('default');
    else
        rng(randomState);
    end
    
    if istable(X)
        X = table2array(X); 
    end
    
    if istable(Y)
        Y = table2array(Y); 
    end
    
    % Shuffle the dataset using the cvpartition function
    cvObj = cvpartition(size(X, 1), 'HoldOut', trainingSize);
    idx = cvObj.test;
    
    % Split into training and test data based on training size specified
    data.Xtest  = X(~idx, :); data.Xtrain = X(idx,  :);
    data.Ytest  = Y(~idx, :); data.Ytrain = Y(idx,  :);
    
end