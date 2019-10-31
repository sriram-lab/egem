%This function runs the PRIME algorithm
%Input:  1. A model with a given media
%        2. Data - Gene expression data that includes three fields: (a) gene_name - a list of genes
%                  as entrez IDs (a column vector); (b) GE - expression values matrix (genes = rows;
%                  samples = columns); (c) GR - measured growth rates (a column vector)
%        3. min_range, max_range and max_bound as calculated in the previous steps
%Output: 1. GR - The predicted growth rate
%        2. stat - Status of returned solution
%        3. LB_All/UB_all - The lower and upper bounds for each model
function[GR,stat,LB_all,UB_all] = PRIME(model,Data,min_range,max_range,max_bound)

%Saving the original model before splitting
orig_model = model;

%Splitting each reaction to its forward and backward direction
[model,rev_map] = SplitRevRxns (model);
model.ub(model.ub==1000) = max_bound;

%Identify the set of growth associated reactions
[vals,rho,model_rxns] = identifyGrowthAssociatedRxns(model,Data);

%Transform the expression values to bounds acording to the nomalization
%range and direction of correlation
rho = rho./abs(rho);
rho(find(isnan(rho))) = 1;
R = repmat(rho',1,size(vals,2));
Bound = R.*vals;
min_val = ceil(abs(min(min(Bound))));
Bound = Bound+min_val;
for i=1:size(Bound,1)
    Bound(i,:) = ((Bound(i,:)-min(Bound(i,:)))./(max(Bound(i,:))-min(Bound(i,:)))*(max_range-min_range))+min_range;
end


LB = model.lb;
UB = model.ub;
%For each cell-line in the data
for i=1:size(vals,2)
    model.lb = LB;
    model.ub = UB;
    model.ub(model_rxns) =Bound(:,i);%Changing the bounds on the selected set of reactions
    
    %reconstructing the original structure of the model to include both
    %reversible and unidirectional reactions
    [model] = unionModel(model,orig_model,max_bound,model_rxns,rev_map);
    
    %The specific model bounds
    UB_all(:,i) = model.ub;
    LB_all(:,i) = model.lb;
    
    %Computing max_biomass
    Res = RunTomlabLP(model,0);
    GR(i) = Res.result_opt;
    stat(i) = Res.result_status;
end

%This function reconstruct the original structure of the model of both
%bidirectional and unidirectional reactions
function[model] = unionModel(model,orig_model,max_bound,model_rxns,rev_map)

rxns_num = length(orig_model.rxns);
UB = model.ub(model_rxns);
[model_rxns,idx] = sort(model_rxns);
UB = UB(idx);

model = orig_model;

model.ub(model.ub==1000) = max_bound;
model.lb(model.lb==-1000) = -max_bound;

[u,ia,ib] = intersect(model_rxns,rev_map(:,1));
model.lb(rev_map(ib,2)) = -UB(ia);
idx = find(model_rxns>rxns_num);
model.ub(model_rxns(1:idx(1)-1)) = UB(1:idx(1)-1);


