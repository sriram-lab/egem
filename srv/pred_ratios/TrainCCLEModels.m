%% Predict cancer metabolism from histone markers
% *Author*: Scott Campit

%% Summary
% This notebook trains regressors on global chromatin profiles from the Cancer 
% Cell Line Encylopedia to predict metabolite and metabolite ratio levels. Other 
% validation sets are also used from the trained CCLE models.

load /nfs/turbo/umms-csriram/scampit/Data/Models/Ratios/Ratio_data.mat
                  
%% Predicting metabolism from GCP with hyperparameter optimization and K-fold cross validation
% Next, let's perform some hyperparameter optimzation for the regressors using 
% the training set. Note that |regressorEnsemble| is a custom function that performs 
% hyperparameter optimization for the following models:

%% 
% # Robust regression* - if the data supports ROLS
% # Ridge regression
% # LASSO
% # Decision trees
% # Least-squares boosting
% # Random forest
%% 
% Hyperparameter optimization was performed using the |regressorEnsemble.mlx| 
% function, which trains several regression models using k-fold cross validation. 
% The code for this script can be found <https://github.com/ScottCampit/ml/blob/master/MATLAB/src/Ensemble/regressorEnsemble.mlx 
% here>.
% 
% Note that the |regressorEnsemble.mlx| function automatically saves the intermediate 
% models, as well as the final models selected from k-fold cross validation. 

%Set the training size and the random number generator for the
%trainTestSplit function.
%parpool;
randomState  = 'default';
trainingSize = 0.8;
data         = trainTestSplit(GCP_norm, MET_norm, ...
                             trainingSize, ...
                             randomState);

%Train several regressor models to predict individual metabolites
kfold = 3;
filename = 'GCP2Met_Ratio_Models.mat';
save(filename);
k_iter = 1;
col_iter = 1;

GCP2MET_mdls = regressorEnsemble(data.Xtrain, data.Ytrain, ...
                                kfold, ...
                                filename, ...
                                k_iter, col_iter);

%% Predicting GCP from metabolism with hyperparameter optimization and K-fold cross validation
% This next section runs two sets of models: one set of models aims to predict 
% [metabolomics <- epigenomics], while the other set of modesl aims to predict 
% [epigenomics <- metabolomics].

% Set the training size and the random number generator for the
% trainTestSplit function.
data = trainTestSplit(MET_norm, GCP_norm, ...
                      trainingSize, ...
                      randomState);

randomState = 'default';
filename = 'MET2GCP_Ratio_Models.mat';
save(filename);
k_iter = 1;
col_iter = 1;
MET2GCP_mdls = regressorEnsemble(data.Xtrain, data.Ytrain, kfold, filename, k_iter, col_iter);

