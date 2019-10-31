%This function identifies the set of significantly correlated reaction
%expression to growth rate
%Input:  1. A model with a given media
%        2. Data - Gene expression data that includes three fields: (a) gene_name - a list of genes
%                  as entrez IDs (a column vector); (b) GE - expression values matrix (genes = rows;
%                  samples = columns); (c) GR - measured growth rates (a column vector)
%Output: 1. vals - Expression value per reaction across all samples
%        2. rho -  The correlation value of each reaction to the growth rate
%        3. model_rxns - The set of reactions associated with growth     
function[vals,rho,model_rxns] = identifyGrowthAssociatedRxns(model,Data)

%For each reaction, the mean over its catalyzing enzymes (when available) is
%calculated.
for i=1:length(model.rxns)
    genes = find(model.rxnGeneMat(i,:));
    genes_unique = model.genes_unique_map(genes);
    genes_unique_entrez = model.genes_unique(genes_unique);
    [u,ia,ib] = intersect(Data.gene_name,genes_unique_entrez);
    if ~isempty(u) %If the reaction is assoziated with gene found in the data
        if length(u) > 1 %If it is associated with more than one gene then we take the mean
            vals(i,:) = mean(Data.GE(ia,:));
        else
            vals(i,:) = Data.GE(ia,:);%Alternatively, the expression of the single gene is taken
        end
    end
end

%Updating reaction indices and vals to those for which we have information
model_rxns = [1:length(model.rxns)];
empty_idx = find(vals(:,1)==0);
vals(empty_idx,:) = [];
model_rxns(empty_idx) = [];

%Computing the correlation between a reaction expression and measured growth rate
for i=1:size(vals,1)
    [rho(i),p(i)] = corr(vals(i,:)',Data.GR,'type','Spearman');
end

%Correcting for multiple hypothesis using FDR and significance level of
alpha = 0.05;
[pID,pN] = FDR(p,alpha);
if isempty(pID)
	fprintf('There are no reactions that are significantly correlated to the growth rate after correcting for muliple hypothesis using FDR with alpha = %d',alpha);
	return;
 end

%Identifying the set of significantly correlated reactions
[p_sort,idx_sort] = sort(p);
sig_p = find(p_sort<=pID);
idx = idx_sort(1:length(sig_p));

%Update all variables to include only the set found above
vals = vals(idx,:);
rho = rho(idx);
model_rxns = model_rxns(idx);


function [pID,pN] = FDR(p,q)
% FORMAT pt = FDR(p,q)
%
 % p   - vector of p-values
 % q   - False Discovery Rate level
 %
 % pID - p-value threshold based on independence or positive dependence
 % pN  - Nonparametric p-value threshold
 %______________________________________________________________________________
 % @(#)FDR.m    1.3 Tom Nichols 02/01/18
 
 p = sort(p(:));
 V = length(p);
 I = (1:V)';
 
 
 cVID = 1;
 cVN = sum(1./(1:V));
 
 pID = p(max(find(p<=I/V*q/cVID)));
 pN = p(max(find(p<=I/V*q/cVN)));