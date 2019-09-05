initCobraToolbox
changeCobraSolver('gurobi');

%% Tissue Analysis
tissueNoCompSoln = tissueAnalysis('NoComp');
tissueCompSoln = tissueAnalysis('Comp');
tissueFVASoln = tissueAnalysis('FVA');

%% Media Analysis
mediaNoCompSoln = mediaAnalysis('NoComp');
mediaCompSoln = mediaAnalysis('Comp');
mediaFVASoln = mediaAnalysis('FVA');

%% Culture Analysis
cultureNoCompSoln = cultureAnalysis('NoComp');
cultureCompSoln = cultureAnalysis('Comp');
cultureFVASoln = cultureAnalysis('FVA');