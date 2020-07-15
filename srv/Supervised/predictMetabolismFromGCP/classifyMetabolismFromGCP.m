%% Classifying cancer metabolism from histone markers
%% Summary
% This module aims to predict cancer metabolism targets (ie metabolomics) using 
% the global chromatin profiles from the same dataset (CCLE).
%% Load data
% The |predMetabolismFromGCP.mat| file contains the metabolomics (Y) and GCPs 
% (X) mapped with intersecting cancer cell lines between both datasets.

clear all; load predMetabolismFromGCP.mat;

p=normpdf(Y);
Y_class = Y; Y_class(p > 0.05) = 0;
Y_class(p < 0.05 & Y_class > 0) = 1;
Y_class(p < 0.05 & Y_class < 0) = -1;
%% Train Regressors to Predict Individual Metabolite Concentrations
% This module trains several regressor models to predict individual metabolomics 
% values based on the global chromatin profiles.
% Split data into a training and validation set
% First, let's split the data into a training and validation set.

% Set the training size and the random number generator for the
% trainTestSplit function.
trainingSize = 0.7;
randomState = 'default';
data = trainTestSplit(X, Y_class, ...
                      trainingSize, ...
                      randomState);
% Hyperparameter optimization
% Next, let's perform some hyperparameteResor optimzation using the training 
% set.

% Train several regressor models to predict individual metabolites
mdl = classifierEnsemble(data.Xtrain, data.Ytrain);
save('trainedClassifiersAll.mat')
% Evaluating model performances on regression for individual metabolites
% Finally, we'll evaluate the models using metrics such as correlation and error.

load trainedClassifiers.mat

% Evaluate several model metrics 
axis = 1;
metrics_metab = evaluateRegressors(mdl, data, axis);
% Evaluating model performances on regression for individual cancer cell lines
% Finally, we'll evaluate the models using metrics such as correlation and error.

% Evaluate several model metrics 
axis = 0;
metrics_ccl = evaluateRegressors(mdl, data, axis);
%% Visualize the metrics
% Metabolite predictions
% Get dataset ready for visualizations

for i = 1:length(metrics_metab.pearson)
    metab_pearson(i, :) = metrics_metab.pearson{i};
    metab_ppval(i, :)   = metrics_metab.pearspval{i};
    metab_r2(i, :)      = metrics_metab.pearsR2{i};
    metab_mae(i, :)     = metrics_metab.meanAbsSqr{i};
    metab_mse(i, :)     = metrics_metab.meanSqrErr{i};
end
metab_pearson = metab_pearson';
metab_mae = metab_mae';
metab_mse = metab_mse';

%%
% Boxplots
%figure;
%boxplot(ccl_pearson, 'orientation', 'horizontal'); xline(0);
% Cancer cell line predictions
% Get dataset ready for visualizations

% Get a single array for each metric
for i = 1:length(metrics_ccl.pearson)
    ccl_pearson(i, :) = metrics_ccl.pearson{i};
    ccl_ppval(i, :)   = metrics_ccl.pearspval{i};
    ccl_r2(i, :)      = metrics_ccl.pearsR2{i};
    ccl_mae(i, :)     = metrics_ccl.meanAbsSqr{i};
    ccl_mse(i, :)     = metrics_ccl.meanSqrErr{i};
end
ccl_pearson = ccl_pearson';
ccl_mae = ccl_mae';
ccl_mse = ccl_mse';
save('trainedClassifiersAllMetrics.mat')
% Heatmap showing correlations and errors 

% Hide non-sigificant correlations
%mask = ccl_ppval > 0.05;
%tmp = ccl_pearson;
%tmp(mask) = NaN;
% Visualize in different ways
%modelNames = ["Robust LS", "Ridge", "LASSO", "DT", "Boosting", "RF"];
%heatmap(tmp, 'YData', modelNames)
%% Use histone markers that have a known metabolic effect
% This module trains several regressor models to predict individual metabolomics 
% values based on the global chromatin profiles that are known to affect metabolism.
% Split data into a training and validation set
% First, let's split the data into a training and validation set.

p=normpdf(Y_metab);
Y_class = Y; Y_class(p > 0.05) = 0;
Y_class(p < 0.05 & Y_class > 0) = 1;
Y_class(p < 0.05 & Y_class < 0) = -1;

% Set the training size and the random number generator for the
% trainTestSplit function.
trainingSize = 0.7;
randomState = 'default';
data = trainTestSplit(X_metab, Y_class, ...
                      trainingSize, ...
                      randomState);
% Hyperparameter optimization
% Next, let's perform some hyperparameter optimzation using the training set.

% Train several regressor models to predict individual metabolites
metabMdl = ensemble(data.Xtrain, data.Ytrain);
save('trainedClassifiersMetabOnly.mat')
% Evaluating model performances on regression for individual metabolites
% Finally, we'll evaluate the models using metrics such as correlation and error.

load trainedClassifiersMetabOnly.mat

% Evaluate several model metrics 
axis = 1;
metrics_metab = evaluateRegressors(mdl, data, axis);
% Evaluating model performances on regression for individual cancer cell lines
% Finally, we'll evaluate the models using metrics such as correlation and error.

% Evaluate several model metrics 
axis = 0;
metrics_ccl = evaluateRegressors(mdl, data, axis);
%% Visualize the metrics
% Metabolite predictions
% Get dataset ready for visualizations

for i = 1:length(metrics_metab.pearson)
    metab_pearson(i, :) = metrics_metab.pearson{i};
    metab_ppval(i, :)   = metrics_metab.pearspval{i};
    metab_r2(i, :)      = metrics_metab.pearsR2{i};
    metab_mae(i, :)     = metrics_metab.meanAbsSqr{i};
    metab_mse(i, :)     = metrics_metab.meanSqrErr{i};
end
metab_pearson = metab_pearson';
metab_mae = metab_mae';
metab_mse = metab_mse';

%%
% Boxplots
%figure;
%boxplot(ccl_pearson, 'orientation', 'horizontal'); xline(0);
% Cancer cell line predictions
% Get dataset ready for visualizations

% Get a single array for each metric
for i = 1:length(metrics_ccl.pearson)
    ccl_pearson(i, :) = metrics_ccl.pearson{i};
    ccl_ppval(i, :)   = metrics_ccl.pearspval{i};
    ccl_r2(i, :)      = metrics_ccl.pearsR2{i};
    ccl_mae(i, :)     = metrics_ccl.meanAbsSqr{i};
    ccl_mse(i, :)     = metrics_ccl.meanSqrErr{i};
end
ccl_pearson = ccl_pearson';
ccl_mae = ccl_mae';
ccl_mse = ccl_mse';
% Heatmap showing correlations and errors 
save('trainedClassifiersMetabOnlyMetrics.mat')
% Hide non-sigificant correlations
%mask = ccl_ppval > 0.05;
%tmp = ccl_pearson;
%tmp(mask) = NaN;
% Visualize in different ways
%modelNames = ["Robust LS", "Ridge", "LASSO", "DT", "Boosting", "RF"];
%heatmap(tmp, 'YData', modelNames)
