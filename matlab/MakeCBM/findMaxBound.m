%The function identifies a global new maximal bound that do not effect cell
%growth. The procedure is done by gradually reducing the bounds on all reactions at once
%
%Input - A model with a given media
%Output - maximal bound
function[max_bound] = findMaxBound(model)

%Splitting each reaction to its forward and backward direction
[model] = SplitRevRxns (model);
Res = RunTomlabLP(model,1);
max_biomass = Res.result_opt;

eps = 1e-4;
idx = find(model.ub==1000);
vals = [1000:-1:0];
%Gradually reducing the upper bound of all reactions in steps of 0.1.
for i=1:length(vals)
    model.ub(idx) = vals(i);
    Res = RunTomlabLP(model,0);
    max_b = Res.result_opt;
    if max_b<max_biomass-eps
        max_bound = vals(i);
        break;
    end
end