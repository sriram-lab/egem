function model = regressorEnsemble(X, Y, kfold, filename, k_iter, col_iter)
%% REGRESSORENSEMBLE Training an ensemble of regressors
% *Author*: Scott Campit
% 
% This function runs an ensemble of regressors for prototyping, although these 
% models are fine-tuned using k-fold cross validation and Bayesian hyperparameter 
% optimization. 
% 
% Currently, there are 6 statistical learning models that are encoded in this 
% function:
%% 
% # Robust least squares
% # Ridge regression
% # LASSO
% # Decision trees
% # Boosting
% # Random forests
%% 
% Note that a dependency in this script is MATLAB's Parallel Computing Toolbox. 
% If it is not possible to download, comment out the |'UseParallel', true, ...| 
% argument in the script.  
% 
% *INPUTS*=
% 
% |X|: A numerical array containing the features as columns and samples as rows.
% 
% |Y|: A numerical array containing the targets as columns and samples as rows.
% 
% |kfold (optional)|: A scalar value denoting the number of k-splits for the 
% dataset.
% 
% *OUTPUT*
% 
% |model|: A structure of several machine learning models from hyperparameter 
% optimization.
    
    try % Convert to GPU array if possible. This may speed up computations downstream.
        X = gpuarray(X);
        Y = gpuarray(Y);
    catch ME
    end
    
    if isempty(kfold)
        kfold = 5;
    end
    
    if isempty(k_iter)
        k_iter = 1;
    end
    
    if isempty(col_iter)
        col_iter = 1;
    end
    
    
    diary regressor.log
    
    data2 = trainTestSplit(X, Y, 0.8, 'Default');
    Xtrain2 = data2.Xtrain;
    Ytrain2 = data2.Ytrain;
    Xtest2  = data2.Xtest;
    Ytest2  = data2.Ytest;
        
    % Cross validation
    for i = k_iter:kfold
        
        % Generate cross validation indices from the Ytrain variable.
        idx = crossvalind('Kfold', size(Y, 1), kfold);        
        
        % Create cross-validation sets
        Xtrain2 = X(idx ~= i, :);
        Ytrain2 = Y(idx ~= i, :);
        Xtest2  = X(idx == i, :);
        Ytest2  = Y(idx == i, :);
        % Train and evaluate univariate models
        for j = col_iter:size(Y, 2) % Uses the parallel toolbox for hyperparameter optimization

                 % Robust Least Squares (If Possible)
                 tmp_lm                = fitlm(Xtrain2, Ytrain2(:, j), ...
                                                'RobustOpts','on'); 
                 coefs                  = table2array(tmp_lm.Coefficients(2:end, 'Estimate'));
                 ypred2                 = predict(tmp_lm, Xtest2);
                 tmp_acc.lm_mdl(i, j)   = corr(ypred2, Ytest2(:, j));
 
                 if i == 1
                     if ~isempty(coefs)
                        model.lm_mdl(:, j)  = coefs;
                     else
                         model.lm_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
                     end
                 else
                     if abs(tmp_acc.lm_mdl(i, j)) > abs(tmp_acc.lm_mdl(i-1, j))
                         if ~isempty(coefs)
                            model.lm_mdl(:, j) = coefs;
                         else
                            model.lm_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
                         end
                     end
                 end

             % Ridge regression
             
                 tmp_ridge                   = fitrlinear(Xtrain2, Ytrain2(:, j), ...
                                                'Learner', 'leastsquares', ...
                                                'Regularization', 'ridge', ...
                                                'OptimizeHyperparameters', 'auto', ...
                                                'HyperparameterOptimizationOptions', ...
                                                  struct('AcquisitionFunctionName', ...
                                                         'expected-improvement-plus', ...
                                                         'UseParallel', true, ...
                                                         'Verbose', 0, ...
                                                         'ShowPlots', false));
                 coefs                      = tmp_ridge.Beta;
                 ypred2                     = predict(tmp_ridge, Xtest2);
                 tmp_acc.ridge_mdl(i, j)    = corr(ypred2, Ytest2(:, j));
 
                 if i == 1
                     if ~isempty(coefs)
                        model.ridge(:, j)  = coefs;
                     else
                         model.ridge(:, j) = zeros([size(Xtrain2, 2), 1]);
                     end
                 else
                     if abs(tmp_acc.ridge_mdl(i, j)) > abs(tmp_acc.ridge_mdl(i-1, j))
                         if ~isempty(coefs)
                            model.ridge(:, j)  = coefs;
                         else
                            model.ridge(:, j) = zeros([size(Xtrain2, 2), 1]);
                         end
                     end
                 end

             
             % LASSO

		 tmp_lasso                     = fitrlinear(Xtrain2, Ytrain2(:, j), ...
		                                'Learner', 'leastsquares', ...
		                                'Regularization', 'lasso', ...
		                                'OptimizeHyperparameters', 'auto', ...
		                                'HyperparameterOptimizationOptions', ...
		                                  struct('AcquisitionFunctionName', ...
		                                         'expected-improvement-plus', ...
		                                         'UseParallel', true, ...
		                                         'Verbose', 0, ...
		                                         'ShowPlots', false));
		 coefs                      = tmp_lasso.Beta;
		 ypred2                     = predict(tmp_lasso, Xtest2);
		 tmp_acc.lasso_mdl(:, j)    = corr(ypred2, Ytest2(:, j));

		 if i == 1
             if ~isempty(coefs)
                model.lasso_mdl(:, j)  = coefs;
             else
                 model.lasso_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
             end
		 else
		     if (abs(tmp_acc.lasso_mdl(i, j)) > abs(tmp_acc.lasso_mdl(i-1, j))) | (all(model.lasso_mdl(:, j)))
                 if ~isempty(coefs)
                    model.lasso_mdl(:, j)  = coefs;
                 else
                     model.lasso_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
                 end
		     end
		 end
     
%             % Decision tree regressors                                  
%             tmp_mdl                    = fitrtree(Xtrain2, Ytrain2(:, j), ...
%                                         'OptimizeHyperparameters', 'auto', ...
%                                         'HyperparameterOptimizationOptions', ...
%                                              struct('AcquisitionFunctionName', ...
%                                                     'expected-improvement-plus', ...
%                                                     'UseParallel', true, ...
%                                                     'Verbose', 0, ...
%                                                     'ShowPlots', false));
%             ypred2                     = predict(tmp_mdl, Xtest2);
%             tmp_acc.dtr_mdl(i, j)      = corr(ypred2, Ytest2(:, j));   
%             
%             if i == 1
%                 model.dtr_mdl{j} = tmp_mdl;
%             else
%                 if abs(tmp_acc.dtr_mdl(i, j)) > abs(tmp_acc.dtr_mdl(i-1, j))
%                     model.dtr_mdl{j}  = tmp_mdl;
%                 end
%             end
%             
%             % Least-squares boosting
%             tmp_mdl                    = fitrensemble(Xtrain2, Ytrain2(:, j), ...
%                                              'Method', 'LSBoost', ...
%                                              'OptimizeHyperparameters', {'NumLearningCycles', ...
%                                                                     'LearnRate'}, ...
%                                              'HyperparameterOptimizationOptions', ...
%                                              struct('AcquisitionFunctionName', ...
%                                                     'expected-improvement-plus', ...
%                                                     'UseParallel', true, ...
%                                                     'Verbose', 0, ...
%                                                     'ShowPlots', false));
%             ypred2                     = predict(tmp_mdl, Xtest2);
%             tmp_acc.lsb_mdl(i, j)      = corr(ypred2, Ytest2(:, j));
%              
%             if i == 1
%                 model.lsb_mdl{j}  = tmp_mdl;
%             else
%                 if abs(tmp_acc.lsb_mdl(i, j)) > abs(tmp_acc.lsb_mdl(i-1, j))
%                     model.lsb_mdl{j}  = tmp_mdl;
%                 end
%             end
%             
            % Random forest regressors
%            try
%                tmp_mdl                    = fitrensemble(Xtrain2, Ytrain2(:, j), ...
%                                                'Method', 'Bag', ...
%                                                'OptimizeHyperparameters', {'NumLearningCycles', ...
%                                                                        'MinLeafSize', ...
%                                                                        'MaxNumSplits'}, ...
%                                                 'HyperparameterOptimizationOptions', ...
%                                                 struct('AcquisitionFunctionName', ...
%                                                        'expected-improvement-plus', ...
%                                                        'UseParallel', true, ...
%                                                        'Verbose', 0, ...
%                                                        'ShowPlots', false));
%                ypred2                     = predict(tmp_mdl, Xtest2);
%                tmp_acc.rfr_mdl(i, j)      = corr(ypred2, Ytest2(:, j));
%
%                if i == 1
%                    model.rfr_mdl{j}  = tmp_mdl;
%                else
%                    if abs(tmp_acc.rfr_mdl(i, j)) > abs(tmp_acc.rfr_mdl(i-1, j))
%                        model.rfr_mdl{j}  = tmp_mdl;
%                    end
%                end
%            catch ME
%                mode.rfr_mdl{j} = NaN;
%            end
            
             % Support vector regressors
		 tmp_svm                   = fitrsvm(Xtrain2, Ytrain2(:, j), ...
			                         'OptimizeHyperparameters', 'auto', ...
			                          'HyperparameterOptimizationOptions', ...
			                          struct('AcquisitionFunctionName', ...
			                                 'expected-improvement-plus', ...
			                                 'UseParallel', true, ...
			                                 'Verbose', 0, ...
			                                 'ShowPlots', false));
		 coefs                      = tmp_svm.Beta;
		 ypred2                     = predict(tmp_svm, Xtest2);
		 tmp_acc.svm_mdl(i, j)      = corr(ypred2, Ytest2(:, j));

		 if i == 1
             if ~isempty(coefs)
                model.svm_mdl(:, j)  = coefs;
             else
                 model.svm_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
             end
		 else
		     if (abs(tmp_acc.svm_mdl(i, j)) > abs(tmp_acc.svm_mdl(i-1, j))) | (all(model.svm_mdl(:, j)))
                 if ~isempty(coefs)
                    model.svm_mdl(:, j)  = coefs;
                 else
                     model.svm_mdl(:, j) = zeros([size(Xtrain2, 2), 1]);
                 end
		     end
         end
        save(filename, 'model', 'i', 'j', '-append');
            
        end
    end
end
