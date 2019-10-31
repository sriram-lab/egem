%% Linear regression for CCLE dataset
CCLEComp = readtable('CCLECompModel.csv');
CCLEFVA = readtable('CCLEFVAModel.csv');

% Clean up models for MATLAB
cellLines = CCLEComp(:, 1);
CCLEComp(:, 1) = [];
CCLEFVA(:, 1) = [];

%% 1. All data
CCLECompAll = modeller(CCLEComp, 'all');
CCLEFVAAll = modeller(CCLEFVA, 'all');

%% 2. Gene expression only
CCLECompGene = modeller(CCLEComp, 'geneExpOnly');
CCLEFVAGene = modeller(CCLEFVA, 'geneExpOnly');

%% 3. Flux data only
CCLECompFlux = modeller(CCLEComp, 'fluxOnly');
CCLEFVFlux = modeller(CCLEFVA, 'fluxOnly');

%% 4. Individual histone marks
CCLECompMe1 = modeller(CCLEComp, 'Me1');
CCLEFVAMe1 = modeller(CCLEFVA, 'Me1');

CCLECompMe2 = modeller(CCLEComp, 'Me2');
CCLEFVMe2 = modeller(CCLEFVA, 'Me2');

CCLECompMe3 = modeller(CCLEComp, 'Me3');
CCLEFVAMe3 = modeller(CCLEFVA, 'Me3');

CCLECompAc = modeller(CCLEComp, 'Ac');
CCLEFVAc = modeller(CCLEFVA, 'Ac');

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

