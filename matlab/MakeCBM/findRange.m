%This function identifies the minimal and maximal values of the
%normalization range
%Input: 1. A model with a given media
%       2. A set of essential reactions as indices in the split model (computed over the split model; a column vector)
%       3. The max_bound identified in the previous step
%       4. Data - Gene expression data that includes three fields: (a) gene_name - a list of genes
%                 as entrez IDs (a column vector); (b) GE - expression values matrix (genes = rows, samples = columnns); 
%                (c) GR - measured growth rates (a column vector)
%Output: minimal and maximal range
function[min_range,max_range] = findRange(model,essential_rxns,max_bound,Data)

%Splitting each reaction to its forward and backward direction
[model] = SplitRevRxns (model);

%Setting the new maximal bound
model.ub(model.ub==1000) = max_bound;

%Computing maximal biomass and setting the biomass reaction to 10% of its maximum
Res = RunTomlabLP(model,1);
max_biomass = Res.result_opt;
biomass = find(model.c==1);
model.lb(biomass) = 0.1*max_biomass;

%Calculating the min_range - the maximal of all minimal values necessary
%to produce 10% maximal biomass
for i=1:length(essential_rxns)
    model.c = zeros(length(model.rxns),1);
    model.c(essential_rxns(i)) = -1;
    Res = RunTomlabLP(model,0);
    min_val(i) = Res.result_opt;
end

min_range = max(min_val);

%Re-setting biomass bound
model.lb(biomass) = 0;
model.c = zeros(length(model.rxns),1);
model.c(biomass) = 1;

counter=1;
[vals,rho,model_rxns] = identifyGrowthAssociatedRxns(model,Data);

%Computing the effect of bound change on growth rate
for i=min_range:0.1:max_bound
    model.ub(model_rxns) = i;
    Res = RunTomlabLP(model,0);
    GR(counter) = Res.result_opt;
    counter = counter+1;
end

%Computing the max_range by searching for the highest point in the segment
%that mostly affect growth
counter=1;
for i=1:length(GR)-1
    diff(counter) = abs(GR(i)-GR(i+1));
    counter = counter+1;
end

counter=1;
for i=1:length(diff)-1
    num = abs(diff(i)-diff(i+1));
    diff2(counter) = str2num(sprintf('%3.5f',num));%Dealing with some sensitivity issues.
    counter = counter+1;
end

[a,i] = max(diff2);
vals = [min_range:0.1:max_bound];
max_range = vals(i+1);


