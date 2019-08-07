function test_media()

load ./../models/recon1
model = metabolicmodel;

queried_medium = 'RPMI';
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

disp("Success!")

end

