% @author: Scott Campit
function [TableChecks] = testMetabolicModel(model)
    changeCobraSolver('gurobi');
    testModel = model;
    
    EXChecks = testEXLeaks(testModel);
    DMChecks = testDMLeaks(testModel);
    DuplicateChecks = testDuplicateReactions(testModel);
    emptyRxnColumnPosition = testRxnGeneMat(testModel);
    wrongLB = testDMRxnLB(testModel);
    
    % TEST 7: Identify Dead End Metabolites
    outputMets = detectDeadEnds(model)
    DeadEnds = model.mets(outputMets)
    [rxnList, rxnFormulaList] = findRxnsFromMets(model, DeadEnds)
    
    model.lb(find(ismember(model.rxns,rxnList)))
    model.ub(find(ismember(model.rxns,rxnList)))

    [allGaps, rootGaps, downstreamGaps] = gapFind(model, 'true')
    
    % TEST 8: Idenitfy Blocked Reactions

    end