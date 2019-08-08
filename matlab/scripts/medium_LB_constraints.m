function model = medium_LB_constraints(model, queried_medium)
% medium_LB_constraints.m sets the lower bound constraints on the Genome-scale metabolic model. 
% By default, the substrate uptake rates were set to RPMI conditions by default. 
% Other medium conditions were scaled w.r.t RPMI amounts (using ratios
% from concentrations as the scaling factor).
% @author: Scott Campit
    
    path = './../../data/final_medium_conditions.xlsx';
    if verLessThan('matlab', '9.6.0.1072779')
        [~, sheetNames] = xlsfinfo(path);
        for sheets = 1:length(sheetNames)
            if ismember(sheetNames(sheets), queried_medium)
                [adjustedLB, rxn_ids] = xlsread(path, sheetNames{sheets});
                rxn_ids(1,:) = [];
                rxn_ids(:,1) = [];

                for rxn=1:length(rxn_ids)
                    model.lb(find(ismember(model.rxns, rxn_ids(rxn, 1)))) = ...
                        adjustedLB(rxn, 4);
                end

            elseif ismember({'nan'}, queried_medium)
                [adjustedLB, rxn_ids] = xlsread(path, 'RPMI');
                rxn_ids(1,:) = [];
                rxn_ids(:,1) = [];

                for rxn=1:length(rxn_ids)
                    model.lb(find(ismember(model.rxns, rxn_ids(rxn, 1)))) = ...
                        adjustedLB(rxn, 4);
                end

            end
        end

    else
        [~, sheetNames] = xlsfinfo(path);
        for sheets = 1:length(sheetNames)
            if ismember(sheetNames(sheets), queried_medium)
                dataArray = readcell(path,...
                    'Sheet', sheetNames(sheets));
                dataArray(1,:) = [];
                dataArray(:,1) = [];

                for rxn=1:length(dataArray)
                    model.lb(find(ismember(model.rxns, dataArray(rxn, 1)))) = ...
                        cell2mat(dataArray(rxn, 7));
                end

            elseif ismember({'nan'}, queried_medium)
                dataArray = readcell('./../../data/final_medium_conditions.xlsx',...
                    'Sheet', 'RPMI');
                dataArray(1,:) = [];
                dataArray(:,1) = [];

                for rxn=1:length(dataArray)
                    model.lb(find(ismember(model.rxns, dataArray(rxn, 1)))) = ...
                        cell2mat(dataArray(rxn, 7));
                end

            end
        end
    end

end