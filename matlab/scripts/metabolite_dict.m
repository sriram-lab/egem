function [excess, depletion] = metabolite_dict(struct, metabolites, media,...
    sheetname, flag)

% Load substrate uptake rates, medium components, reactions of interest
load ./../vars/supplementary_software_code media_exchange1
var = {...
    './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat' ...
    };

for kk = 1:numel(var)
    load(var{kk});
end

% Need labels for x and y axis 
media = string(mediareactions1(:,2));
rxns = string(metabolites(:,3));

switch flag
    case {'sra', 'competition', 'no_competition'}
        ezscr = normalize(struct.excess_flux);
        ekeep = ezscr > 2;
        dzscr = normalize(struct.depletion_flux);
        dkeep = dzscr > 2;
        [array_row, array_col] = size(struct.excess_flux);

    case 'fva'
        % Create logical with excess and depletion arrays
        ezscr = normalize(struct.excess_maxflux);
        ekeep = ezscr > 2;
        dzscr = normalize(struct.depletion_maxflux);
        dkeep = dzscr > 2;
        [array_row, array_col] = size(struct.excess_maxflux);
end

% Capture medium components that have zscr > 2 w.r.t. reaction (excess)
for i=1:length(media)
    for j = 1:length(rxns)
        if ekeep(i,j) == 1
            excess(i,j) = media(i);
        elseif ekeep(i,j) == 0
            excess(i,j) = "NaN";
        end
        disp(j)
    end
    disp(i)
end

% Capture medium components that have zscr > 2 w.r.t. reaction (depletion)
for i=1:length(media)
    for j = 1:length(rxns)
        if dkeep(i,j) == 1
            depletion(i,j) = media(i);
        elseif dkeep(i,j) == 0
            depletion(i,j) = "NaN";
        end
        disp(j)
    end
    disp(i)
end

% Output dictionary values and save to Excel File
rows = {'Excess', 'Depletion'};
columns = metabolites(:,3);

excess = num2cell(cellstr(excess)',2);
depletion = num2cell(cellstr(depletion)',2);
excess = [excess{:}];
excess = reshape(excess, [array_row, array_col]);
depletion = [depletion{:}];
depletion = reshape(depletion, [array_row, array_col]);

switch flag
    case 'competition' 
        switch media
            case 'RPMI'
                disp(media)
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A1');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A2:A3');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B1:U1');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B2:U2');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B3:U3');
            case 'DMEM'
                disp(media)
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A5');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A6:A7');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B5:U5');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B6:U6');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B7:U7');
            case 'L15'
                disp(media)
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A9');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A10:A11');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B9:U9');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B10:U10');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B11:U11');
        end
    case 'no_competition'
        switch media
            case 'RPMI'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'W1');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'W2:W3');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'X1:AQ1');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'X2:AQ2');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'X3:AQ3');
            case 'DMEM'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'W5');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'W6:W7');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'X5:AQ5');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'X6:AQ6');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'X7:AQ7');
            case 'L15'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'W9');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'W10:W11');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'X9:AQ9');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'X10:AQ10');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'X11:AQ11');
        end
    case 'fva'
        switch media
            case 'RPMI'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A1');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A2:A3');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B1:U1');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B2:U2');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B3:U3');
            case 'DMEM'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A5');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A6:A7');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B5:U5');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B6:U6');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B7:U7');
            case 'L15'
                xlswrite('./metabolic_sensitivity.xlsx', strcat(flag,'_',media), sheetname, 'A9');
                xlswrite('./metabolic_sensitivity.xlsx', rows, sheetname, 'A10:A11');
                xlswrite('./metabolic_sensitivity.xlsx', columns, sheetname, 'B9:U9');
                xlswrite('./metabolic_sensitivity.xlsx', excess, sheetname, 'B10:U10');
                xlswrite('./metabolic_sensitivity.xlsx', depletion, sheetname, 'B11:U11');
        end
end

% % Output dictionary object and save to json file
% excess = num2cell(cellstr(excess)',2);
% depletion = num2cell(cellstr(depletion)',2);
% excess = containers.Map(cellstr(rxns), excess, 'UniformValues', false);
% depletion = containers.Map(cellstr(rxns), depletion, 'UniformValues', false);
% save_json(jsonencode(excess), filename,'_',flag,'_excess_.json'))
% save_json(jsonencode(depletion), filename,'_',flag,'_depletion_.json'))

end
