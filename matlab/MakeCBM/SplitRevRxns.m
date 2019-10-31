%Split all reversible reactions to their forward and backward direction
function [model,rev_map] = SplitRevRxns (model)
rev_rxns = find(model.lb<0);
[num_mets, num_rxns] = size(model.S);
model.rev_map(:,1) = [num_rxns+1:num_rxns+length(rev_rxns)];
model.rev_map(:,2) = rev_rxns;
rev_map = model.rev_map;
model.S = [model.S -model.S(:,rev_rxns)];
orig_lb = model.lb (rev_rxns);
model.lb (rev_rxns) = 0;
model.rev (rev_rxns) = 0;
model.rev = [model.rev; zeros(length(rev_rxns),1)];
model.lb = [model.lb; zeros(length(rev_rxns),1)];
model.ub = [model.ub; -orig_lb];
bkwd_rxns = model.rxns(rev_rxns);
bkwd_rxnNames = model.rxns(rev_rxns);
for i=1:length(rev_rxns)
        name = model.rxns{rev_rxns(i)};
        name_long = model.rxnNames{rev_rxns(i)};
        model.rxns{rev_rxns(i)} = sprintf ('%s_fwd',name);
        bkwd_rxns{i} = sprintf ('%s_bkwd',name);
        model.rxnNames{rev_rxns(i)} = sprintf ('%s_fwd',name_long);
        bkwd_rxnNames{i} = sprintf ('%s_bkwd',name_long);
end
model.rxns = [model.rxns; bkwd_rxns];
model.rxnNames = [model.rxnNames; bkwd_rxnNames];
model.rxnGeneMat = [model.rxnGeneMat; model.rxnGeneMat(rev_rxns,:)];
model.rules = [model.rules; model.rules(rev_rxns)];
model.grRules = [model.grRules; model.grRules(rev_rxns)];
model.subSystems = [model.subSystems; model.subSystems(rev_rxns)];
model.c = [model.c; model.c(rev_rxns)];


end
