% Make eGEM adds reactions for epigenome-scale metabolic modeling
% @author: Scott Campit & Lauren Fane
function [eGEM] = make_eGEM(metabolicmodel)

    %if (~exist('metabolicModel','var')) || (isempty(metabolicmodel))
        % The AcGEM from Shen et al., 2019 is the default base metabolic model. 
        % It was edited by @scampit with the following changes:
        % An extra exchange reaction that is not in the BiGG RECON1 metabolic 
            % was removed. The biomass objective function and the {'DM_Tyrggn_c'}
            % reaction were kept.
        % The direct exchange reaction for peplys_e -> peplys_n was removed
        % Additional database identifiers and metabolite charges were added to 
            % the metabolic model.
    %    load './../../metabolic_models/recon1.mat'
    temporary_model = metabolicmodel;
    %end

    %reaction_path = '/home/scampit/Desktop/eGEM/data/metabolicModel_maps/metabolic_map.xlsx';
    reaction_path = './../../../data/metabolicModel_maps/metabolic_map.xlsx';
    [RxnVals, RxnTxt] = xlsread(reaction_path, 'Reactions');
    reactionIDs = RxnTxt(2:end, 2);
    reactionFormulas = RxnTxt(2:end, 4);
    gprAssociations = RxnTxt(2:end, 9);
    rxnUB = RxnVals(:, 6);
    rxnLB = RxnVals(:, 7);
    rxnRev = RxnVals(:, 5);
    
    for i = 1:length(RxnVals(:,1))
        temporary_model = addReaction(temporary_model, char(reactionIDs(i,1)), ...
            'reactionFormula', char(reactionFormulas(i,1)), ...
            'lowerBound', rxnLB(i,1), ...
            'upperBound', rxnUB(i,1), ...
            'geneRule', char(gprAssociations(i,1)), ...
            'reversible', rxnRev(i,1) ... % Doesn't work with num -> do manually
            );
    end

    [MetVals, MetTxt] = xlsread(reaction_path, 'Metabolites');
    metNames = MetTxt(2:end, 2);
    metFormula = MetTxt(2:end, 3);
    keggID = MetTxt(2:end, 6);
    hmdbID = MetTxt(2:end, 7);
    chebiID = MetTxt(2:end, 10);
    metCharge = MetVals;
    
    original_metLength = length(metabolicmodel.mets);

    for j = 1:length(MetVals)
        temporary_model.metCharge(original_metLength+j) = metCharge(j);
        temporary_model.metNames(original_metLength+j) = metNames(j);
        temporary_model.metFormulas(original_metLength+j) = metFormula(j);
        temporary_model.metChEBIID(original_metLength+j) = chebiID(j);
        temporary_model.metKEGGID(original_metLength+j) = keggID(j);
        temporary_model.metHMDBID(original_metLength+j) = hmdbID(j);
    end
    
    eGEM = temporary_model;
    save('./../../metabolic_models/eGEM.mat', 'eGEM');

    %clear all;
end
