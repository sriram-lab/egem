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
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end
% If medium condition is not specified, set to RPMI
elseif ismember({'nan'}, medium)
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

%% DMEM
    % https://www.thermofisher.com/order/catalog/product/12491015
    % No biotin or vitamin B12
    % Contains 2.25x more D-glc than RPMI
    % Also contains excess of basically everything
elseif ismember({'DMEM'}, medium)
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end
    

%% L15
    % https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.80.html
    % Supports monkey kidney cells (HEP-2) and primary explants of
    % embryonic and adult human tissue
    % includes galactose, phenol red, glutamine and sodium pyruvate
elseif ismember({'L15'}, medium)
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'L15');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.83.html
elseif ismember({'McCoy 5A'} , medium) % McCoy 5A
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'McCoy 5A');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.76.html
elseif ismember({'Iscove'} , medium) % IMDM
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'Iscove');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.126.html
elseif ismember({'Waymouth'}, medium) % Waymouth
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'Waymouth');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.227.html
elseif ismember({'DMEM:F12 (1:1)'}, medium) % 1:1 DMEM and F12
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM-F12');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.64.html
elseif ismember({'HAM F-12'}, medium) % F12
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'HAM F-12');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/order/catalog/product/12571071?SID=srch-srp-12571071
% It is different from DMEM, but need to encode glutamine and nucleosides
% into the code...
elseif ismember({'alpha-MEM'}, medium) % Alternative form to MEM
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'alpha-MEM');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
elseif ismember({'RPMI w Gln'}, medium) % RPMI ++ L-gln
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI w Gln');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.61.html
elseif ismember({'HAM F-10'}, medium) % F10
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'HAM F10');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% No one sells it - calc'd ratio
elseif ismember({'DMEM:RPMI (2:1)'}, medium) % 2:1 DMEM:RPMI
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM-RPMI 2-1');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% No one sells it - calc'd ratio
elseif ismember({'MCDB105:M199'}, medium) % MCDB 105: Medium 199
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'MCDB105-M199');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.314.html
elseif ismember({'Williams E Medium'}, medium) % Williams E Medium
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'Williams E');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% Same as DMEM:F12 
elseif ismember({'ACL-4'}, medium) % ACL-4
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM-F12');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end  

% No one sells it - calc'd ratio
elseif ismember({'RPMI:F12'}, medium) % 1:1 RPMI and F12
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI-F12');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end   

% No one sells it - calc'd ratio
elseif ismember({'DMEM:Iscove'}, medium) % 1:1 DMEM and Iscove 
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM-Iscove');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end

% No one sells it - calc'd ratio  
elseif ismember({'RPMI:Iscove'}, medium) % 1:1 RPMI and Iscove
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI-Iscove');
    for rxn=2:length(txt)
        model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
    end
end

end