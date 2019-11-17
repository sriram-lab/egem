%% Multivariate linear regressiokn
CCLEComp = readtable('CCLECompModel.csv');
CCLEFVA = readtable('CCLEFVAModel.csv');

% Clean up models for MATLAB
cellLines = CCLEComp(:, 1);
CCLEComp(:, 1) = [];
CCLEFVA(:, 1) = [];

array = table2array(CCLEComp);
array = meanCenterPredictors(array);
X0 = array(:, 1:201);
Y = array(:, 202:end);

% beta is the estimated coefficients
% Sigma is the variance matrix of Y
% E is the residuals
% CovB is the variance matrix of the regression coefficients
% logL is the log likeloihood objective function
[beta, Sigma, E, CovB_old, logL] = mvregress(X0, Y);

[betaHat, Yhat, error, CovResidual] = multiVarRegress(X0, Y);
[Xtrain, Ytrain, Xtest, Ytest] = trainTestSplit(X0, Y, 0.8);
ypred = predict(betaHat, Xtest)
alpha = 0.01;
iterations = 100;
[theta, costs] = gradientDescent(X0, Y, beta, alpha, iterations);
