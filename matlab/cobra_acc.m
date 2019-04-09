x = findRxnFromCompartment(model, '[n]');
xlswrite('reactions.xlsx', x);

y = findGenesFromRxns(model, x(:,1));
y = vertcat(y{:})
for i:length(y)
    xlswrite('genes.xlsx', y{i}, ['']
xlswrite('genes.xlsx',y);