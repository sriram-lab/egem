function [excess, depletion] = metabolite_dict(struct, filename, flag)

% Load substrate uptake rates, medium components, reactions of interest
load ./../vars/supplementary_software_code media_exchange1
var = {'./../vars/metabolites.mat', './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk});
end

% Need labels
media = string(mediareactions1(:,2));
rxns = string(metabolites(:,3));

switch flag
    case 'fba'
        ezscr = normalize(struct.excess_flux);
        ekeep = ezscr > 2;
        dzscr = normalize(struct.depletion_flux);
        dkeep = ezscr > 2;
    case 'fva'
        % Create logical with excess and depletion arrays
        ezscr = normalize(struct.excess_max_flux);
        ekeep = ezscr > 2;
        dzscr = normalize(struct.depletion_max_flux);
        dkeep = dzscr > 2;
end

% Capture medium components that have zscr > 2 w.r.t. reaction
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

% Capture medium components that have zscr > 2 w.r.t. reaction
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

% Output dictionary object and save to json file
excess = num2cell(cellstr(excess)',2);
depletion = num2cell(cellstr(depletion)',2);
excess = containers.Map(cellstr(rxns), excess, 'UniformValues', false);
depletion = containers.Map(cellstr(rxns), depletion, 'UniformValues', false);
save_json(jsonencode(excess), strcat(filename,'_',flag,'_excess_.json'))
save_json(jsonencode(depletion), strcat(filename,'_',flag,'_depletion_.json'))

end
