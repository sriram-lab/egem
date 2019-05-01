histone_corr(model, 'amet', [], 'n', 1, 1E-2, 1, 1E-3, 0);

rxnpos1  = [find(ismember(model.rxns, 'EX_KAC'));];
test = read_json('reactions.json')
data = make_heatmap(model, rxnpos1, 1);
