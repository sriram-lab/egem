function model = medium_LB_constraints(model, queried_medium)
% medium_LB_constraints.m sets the lower bound constraints on the Genome-scale metabolic model. 
% By default, the substrate uptake rates were set to RPMI conditions by default. 
% Other medium conditions were scaled w.r.t RPMI amounts (using ratios
% from concentrations as the scaling factor).
% @author: Scott Campit
    
    %path = '/home/scampit/Desktop/eGEM/data/Medium_Component_Maps/final_medium2.xlsx';
    path = 'C:\Users\scampit\Desktop\eGEM\data\Medium_Component_Maps\final_medium2.xlsx';
    if verLessThan('matlab', '9.6.0.1072779')
        [~, sheetNames] = xlsfinfo(path);
        for sheets = 1:length(sheetNames)
            if ismember(string(sheetNames(sheets)), queried_medium)
                [adjustedLB, rxn_ids] = xlsread(path, string(sheetNames{sheets}));
                rxn_ids(1,:) = [];
                rxn_ids(:,1) = [];

                for rxn=1:length(rxn_ids)
                    model.lb(find(ismember(string(model.rxns), string(rxn_ids(rxn, 5))))) = ...
                        adjustedLB(rxn, 2);
                end

            elseif ismember({'nan'}, queried_medium)
                [adjustedLB, rxn_ids] = xlsread(path, 'RPMI');
                rxn_ids(1,:) = [];
                rxn_ids(:,1) = [];

                for rxn=1:length(rxn_ids)
                    model.lb(find(ismember(string(model.rxns), string(rxn_ids(rxn, 5))))) = ...
                        adjustedLB(rxn, 2);
                end

            end
        end

    else
        [~, sheetNames] = xlsfinfo(path);
        for sheets = 1:length(sheetNames)
            if ismember(string(sheetNames(sheets)), string(queried_medium))
                dataArray = readcell(path, 'Sheet', string(sheetNames(sheets)));
                dataArray(1,:) = [];
                dataArray(:,1) = [];

                for rxn=1:length(dataArray)
                    model.lb(find(ismember(string(model.rxns), string(dataArray(rxn, 4))))) = ...
                        cell2mat(dataArray(rxn, 2));
                end

            elseif ismember({'nan'}, queried_medium)
                dataArray = readcell(path, 'Sheet', 'RPMI');
                dataArray(1,:) = [];
                dataArray(:,1) = [];

                for rxn=1:length(dataArray)
                    model.lb(find(ismember(string(model.rxns), string(dataArray(rxn, 4))))) = ...
                        cell2mat(dataArray(rxn, 2));
                end

            end
        end
    end

end