function metrics = evaluateRegressors(model, data, axis)
%% EVALUATEREGRESSORS Evaluate regressors
% *Author*: Scott Campit
% 
% This module evaluates regressor models by computing the Pearson and Spearman 
% correlation coefficient and associated p-values, the coefficient of determination 
% associated with both metrics, and the mean squared error and mean absolute error 
% from the models.
% 
% *INPUTS*
% 
% |model:| A structure containing the regressors to evaluate.
% 
% |data:| A structure containing the test data sets as |Xtest| and |Ytest| respectively.
% 
% *OUTPUT*
% 
% |metrics|: A structure containing model metrics, including correlations and 
% error.
    modelTypes = fieldnames(model);
    Xtest = data.Xtest; Ytest = data.Ytest;
    
    for i = 1:length(modelTypes)
        featMdls = model.(modelTypes{i});
        switch axis
            case 0
                % Compute all ypred
                for j = 1:length(featMdls)
                    featMdl = featMdls{j};
                    ypred(:, j) = predict(featMdl, Xtest);
                end
                % Get R, R2, MAE, MSE
                for k = 1:size(Ytest, 1)
                    covarmat = cov(Ytest(k, :), ypred(k, :));
                    corrmat  = corrcov(covarmat);
                    pears(k) = corrmat(1, 2); % (1, 2) or (2, 1) will work.
                    pr2(k)   = pears(k) .^ 2;
                    mae(k)   = sum(abs(Ytest(k, :) - ypred(k, :)), 'all') ./ size(ypred, 1);
                    mse(k)   = sum((Ytest(k, :) - ypred(k, :)).^2, 'all') ./ size(ypred, 1);
                end
                % Get p-value associated with R
                for l = 1:length(pears)
                    [~, ppval(l)] = ttest2(pears(l), pears);
                end
                spear = []; spval = []; sr2 = [];
            case 1
                for j = 1:length(featMdls)
                    featMdl = featMdls{j};
                    ypred = predict(featMdl, Xtest);
                    Ytest = Ytest; ypred = ypred;
                    [pears(j), ppval(j)] = corr(Ytest(:, j), ypred, 'Type', 'Pearson');
                    [spear(j), spval(j)] = corr(Ytest(:, j), ypred, 'Type', 'Spearman');
                    pr2(j) = pears(j).^2;
                    sr2(j) = spear(j).^2;
                    mae(j) = sum(abs(Ytest(:, j) - ypred), 'all') ./ length(ypred);
                    mse(j) = sum((Ytest(:, j) -ypred).^2, 'all') ./ length(ypred);
                end
        end
        metrics.pearson{i} = pears;  metrics.pearspval{i} = ppval;
        metrics.spreaman{i} = spear; metrics.spearpval{i} = spval;
        metrics.pearsR2{i} = pr2;    metrics.spearR2{i} = sr2;
        metrics.meanAbsSqr{i} = mae; metrics.meanSqrErr{i} = mse;
    end
end