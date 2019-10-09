%% Linear regression for CCLE dataset
CCLEComp = readtable('CCLECompModel.csv');
CCLEFVA = readtable('CCLEFVAModel.csv');

% Clean up models for MATLAB
cellLines = CCLEComp(:, 1);
CCLEComp(:, 1) = [];
CCLEFVA(:, 1) = [];

array = table2array(CCLEComp);
%% 1. All features
% Competition
array = meanCenterPredictors(array);
X0 = array(:, 1:201);
Y = array(:, 202:end);

[beta, Sigma, E, CovB, logL] = mvregress(X0, Y); % Gene expression is the only relevant feature
OutlierPos = find(CCLECompModel.Diagnostics.CooksDistance > ...
3*mean(CCLECompModel.Diagnostics.CooksDistance));
X_clean = X0;
X_clean(OutlierPos, :) = [];
Y_clean = Y;
Y_clean(OutlierPos, :) = [];

CCLECompModel_clean = fitlm(X_clean, Y_clean); 
anova(CCLECompModel_clean)
figure(1)
plotResiduals(CCLECompModel_clean, 'probability');
figure(2)
plotResiduals(CCLECompModel_clean, 'fitted');
qqplot(X_clean(:,1));

% FVA
array = table2array(CCLEFVA);
X0 = array;
X0 = meanCenterPredictors(X0);
Y = X0(:, 2);
X0(:,2) = [];

CCLEFVAModel = fitlm(X0, Y); % Gene expression is the only relevant feature
OutlierPos = find(CCLEFVAModel.Diagnostics.CooksDistance > ...
3*mean(CCLEFVAModel.Diagnostics.CooksDistance));
X_clean = X0;
X_clean(OutlierPos, :) = [];
Y_clean = Y;
Y_clean(OutlierPos, :) = [];

CCLEFVAModel_clean = fitlm(X_clean, Y_clean); 
anova(CCLEFVAModel_clean)
figure(1)
plotResiduals(CCLEFVAModel_clean, 'probability');
figure(2)
plotResiduals(CCLEFVAModel_clean, 'fitted');
qqplot(X_clean(:,1));

%% 2. Gene expression versus flux model
%cov(array)
%rank(array)
%[U, S, V] = svd(array, 'econ');
X1 = array;
X2 = array;

% Leave out gene expression
X1(:,1) = []; 
X1 = meanCenterPredictors(X1);
Y = X1(:, 2);
X1(:, 2) = [];

% Leave out flux data
X2(:, 3:6) = [];
X2 = meanCenterPredictors(X2);
Y = X2(:, 2);
X2(:, 2) = [];

fluxOnlyModel = fitlm(X1, Y);
geneOnlyModel = fitlm(X2, Y);

anova(fluxOnlyModel) % Gene expression is the only relevant feature
anova(geneOnlyModel)

%% 3. One feature only
OneFeatureRegress(X0, Y);
% Results from ANOVA:
    % All features are signficant
    % Gene expression p-value: 2E-17
    % Ac flux p-value: 2.8E-3
    % Me1 flux p-value: 3.0E-50
    % Me2 flux p-value: 2.4E-230
    % Me3 flux p-value: 3.0E-85

% Perform initial data analysis on all data
model1 = fitlm(X, Y);
plotAllDiagnostics(model1)
anova(model1)
plotDiagnostics(model1, 'cookd')
[~, outliers] = max(model1.Diagnostics.CooksDistance);
plotResiduals(model1, 'fitted')
plotResiduals(model1, 'probability')

%% 4. Leave one feature out
leaveOneOutLM(X0, Y)

% Results from ANOVA:
    % Leaving out gene expression results in pvalue = 0
    % leaving out Ac flux results in all signficant features
    % leaving out Me1 flux results in all significant features except Me2
    % flux
    % leaving out Me2 flux results in all significant features
    % leaving out Me3 flux results in all significant features except Me2
    
% This tells me the metabolic flux for Me2 is not a reliable feature

% Remove large leverage points using Cooks distance
newModel1 = fitlm(X, Y, 'Exclude', outliers);
anova(newModel1)
plotDiagnostics(newModel1, 'cookd');
plotResiduals(newModel1, 'fitted');
plotResiduals(newModel1, 'probability');

% Normalize data 
newModel2 = fitlm(log(X), Y, 'Exclude', outliers);
plotResiduals(newModel2, 'probability');
plotResiduals(newModel2, 'fitted');
B = fitrlinear(X, Y);

% Remove transcriptomics data from CCLE
newModel3 = fitlm

%% 5. Histone marker-specific models

CCLEComp = table2cell(CCLEComp);
Me1Comp = CCLEComp(~cellfun(@isempty, regexp(string(CCLEComp(:,3)), 'me1')), :); 
Me2Comp = CCLEComp(~cellfun(@isempty, regexp(CCLEComp(:,3), 'me2')), :); 
Me3Comp = CCLEComp(~cellfun(@isempty, regexp(CCLEComp(:,3), 'me3')), :); 
AcComp = CCLEComp(~cellfun(@isempty, regexp(CCLEComp(:,3), 'ac')), :); 

Me1Labels = Me1Comp(:,1:3);
Me2Labels = Me2Comp(:,1:3);
Me3Labels = Me3Comp(:,1:3);
AcLabels = AcComp(:,1:3);

Me1Comp(:, 1:3) = [];
Me2Comp(:, 1:3) = [];
Me3Comp(:, 1:3) = [];
AcComp(:, 1:3) = [];

X0 = cell2mat(Me2Comp);
X0 = meanCenterPredictors(X0);
Y = X0(:, 2);
Y = normalize(Y, 'range');
X0(:,2) = [];
Me2CompModel = fitlm(X0, Y);

OutlierPos = find(Me2CompModel.Diagnostics.CooksDistance > ...
3*mean(Me2CompModel.Diagnostics.CooksDistance));
X_clean = X0;
X_clean(OutlierPos, :) = [];
Y_clean = Y;
Y_clean(OutlierPos, :) = [];

Me2CompModel_clean = fitlm(X_clean, Y_clean); 
anova(Me2CompModel_clean)

scatter(X_clean(:, 3), Y_clean)
mask = (X_clean(:,3) == min(X_clean(:, 3)));
Me1 = X_clean(~mask, 3);
Y_Me1 = Y(~mask, 1);
scatter(Me1, Y_Me1)