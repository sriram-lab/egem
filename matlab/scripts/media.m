%% @author: Scott Campit
function model = media(model, medium)
%% media.m defines the medium constraints we will impose on the genome-scale metabolic model. 
% By default, the substrate uptake rates were set to RPMI conditions by default. 
    % Other medium conditions were scaled w.r.t RPMI amounts (using ratios
    % from concentrations as the scaling factor).
    
    
%% Amino acid levels (mM/L) in serum of patients and health subjects:
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471778/
    % Alanine
    % Arginine
    % Aspartate
    % Citrate
    % Glutamate
    % Glycine
    % Methionine
    % Phenylalanine
    % Tyrosine
    % Valine

%% RPMI
    % https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
    % Contains glutathione
    % Contains high concentrations of vitamins
        % Contains Biotin
        % Contains Vitamin B12
        % Contains PABA
        % Contains Inositol in high concentrations
        % Contains Choline in high concentrations
if ismember({'RPMI'}, medium)
    num = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMI');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end
% If medium condition is not specified, set to RPMI
elseif ismember({'nan'}, medium)
    num = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMI');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

%% DMEM
    % https://www.thermofisher.com/order/catalog/product/12491015
    % No biotin or vitamin B12
    % Contains 2.25x more D-glc than RPMI
    % Also contains excess of basically everything
elseif ismember({'DMEM'}, medium)
    num = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEM');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end
    

%% L15
    % https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.80.html
    % Supports monkey kidney cells (HEP-2) and primary explants of
    % embryonic and adult human tissue
    % includes galactose, phenol red, glutamine and sodium pyruvate
elseif ismember({'L15'}, medium)
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'L15');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.83.html
elseif ismember({'McCoy 5A'} , medium) % McCoy 5A
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'McCoy 5A');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.76.html
elseif ismember({'Iscove'} , medium) % IMDM
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'Iscove');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.126.html
elseif ismember({'Waymouth'}, medium) % Waymouth
    num = readcell('./../../data/uptake.xlsx', 'Waymouth');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.227.html
elseif ismember({'DMEM:F12 (1:1)'}, medium) % 1:1 DMEM and F12
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'DMEM-F12');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.64.html
elseif ismember({'HAM F-12'}, medium) % F12
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'HAM F-12');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/order/catalog/product/12571071?SID=srch-srp-12571071
% It is different from DMEM, but need to encode glutamine and nucleosides
% into the code...
elseif ismember({'alpha-MEM'}, medium) % Alternative form to MEM
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'alpha-MEM');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
elseif ismember({'RPMI w Gln'}, medium) % RPMI ++ L-gln
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'RPMI w Gln');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.61.html
elseif ismember({'HAM F-10'}, medium) % F10
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'HAM F10');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% No one sells it - calc'd ratio
elseif ismember({'DMEM:RPMI (2:1)'}, medium) % 2:1 DMEM:RPMI
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'DMEM-RPMI 2-1');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% No one sells it - calc'd ratio
elseif ismember({'MCDB105:M199'}, medium) % MCDB 105: Medium 199
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'MCDB105-M199');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.314.html
elseif ismember({'Williams E Medium'}, medium) % Williams E Medium
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'Williams E');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% Same as DMEM:F12 
elseif ismember({'ACL-4'}, medium) % ACL-4
    num = readcell('./../../data/uptake.xlsx', 'Sheet',  'DMEM-F12');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end  

% No one sells it - calc'd ratio
elseif ismember({'RPMI:F12'}, medium) % 1:1 RPMI and F12
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'RPMI-F12');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end   

% No one sells it - calc'd ratio
elseif ismember({'DMEM:Iscove'}, medium) % 1:1 DMEM and Iscove 
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'DMEM-Iscove');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end

% No one sells it - calc'd ratio  
elseif ismember({'RPMI:Iscove'}, medium) % 1:1 RPMI and Iscove
    num = readcell('./../../data/uptake.xlsx',  'Sheet', 'RPMI-Iscove');
    for rxn=2:length(num)
        model.lb(find(ismember(model.rxns, char(num(rxn, 2))))) = cell2mat(num(rxn, 4));
    end
end

end