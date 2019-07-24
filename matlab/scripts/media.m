%% @author: Scott Campit
function model = media(model, medium)
%% media.m defines the medium constraints we will impose on the genomescale metabolic model. 
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
    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.114.html
    % Contains glutathione
    % Contains high concentrations of vitamins
        % Contains Biotin
        % Contains Vitamin B12
        % Contains PABA
        % Contains Inositol in high concentrations
        % Contains Choline in high concentrations
if verLessThan('matlab', '9.6.0.1072779')
    if ismember({'RPMI'}, medium)
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end
    % If medium condition is not specified, set to RPMI
    elseif ismember({'nan'}, medium)
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    %% DMEM
        % https://www.thermofisher.com/order/catalog/product/12491015
        % No biotin or vitamin B12
        % Contains 2.25x more Dglc than RPMI
        % Also contains excess of basically everything
    elseif ismember({'DMEM'}, medium)
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEM');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end


    %% L15
        % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.80.html
        % Supports monkey kidney cells (HEP2) and primary explants of
        % embryonic and adult human tissue
        % includes galactose, phenol red, glutamine and sodium pyruvate
    elseif ismember({'L15'}, medium)
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'L15');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.83.html
    elseif ismember({'McCoy5A'} , medium) % McCoy 5A
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'McCoy5A');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.76.html
    elseif ismember({'Iscove'} , medium) % IMDM
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'Iscove');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.126.html
    elseif ismember({'Waymouth'}, medium) % Waymouth
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'Waymouth');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.227.html
    elseif ismember({'DMEMF12'}, medium) % 1:1 DMEM and F12
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEMF12');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.64.html
    elseif ismember({'HAMF12'}, medium) % F12
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'HAMF12');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/order/catalog/product/12571071?SID=srchsrp12571071
    % It is different from DMEM, but need to encode glutamine and nucleosides
    % into the code...
    elseif ismember({'alphaMEM'}, medium) % Alternative form to MEM
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'alphaMEM');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.114.html
    elseif ismember({'RPMIwGln'}, medium) % RPMI ++ Lgln
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMIwGln');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.61.html
    elseif ismember({'HAMF10'}, medium) % F10
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'HAMF10');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % No one sells it  calc'd ratio
    elseif ismember({'DMEMRPMI21'}, medium) % 2:1 DMEM:RPMI
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEMRPMI21');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % No one sells it  calc'd ratio
    elseif ismember({'MCDB105M199'}, medium) % MCDB 105: Medium 199
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'MCDB105M199');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.314.html
    elseif ismember({'Williams'}, medium) % Williams E Medium
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'Williams');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end

    % Same as DMEM:F12 
    elseif ismember({'ACL4'}, medium) % ACL4
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEMF12');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end 

    % No one sells it  calc'd ratio
    elseif ismember({'RPMIF12'}, medium) % 1:1 RPMI and F12
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMIF12');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end  

    % No one sells it  calc'd ratio
    elseif ismember({'DMEMIscove'}, medium) % 1:1 DMEM and Iscove 
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'DMEMIscove');
        for rxn=2:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn)))) = num(rxn);
        end

    % No one sells it  calc'd ratio  
    elseif ismember({'RPMIIscove'}, medium) % 1:1 RPMI and Iscove
        [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMIIscove');
        txt(1,:) = [];
        for rxn=1:length(txt)
            model.lb(find(ismember(model.rxns, txt(rxn, 2)))) = num(rxn, 2);
        end
    end
    
else
    if ismember({'RPMI'}, medium)
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMI');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end
    % If medium condition is not specified, set to RPMI
    elseif ismember({'nan'}, medium)
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMI');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    %% DMEM
        % https://www.thermofisher.com/order/catalog/product/12491015
        % No biotin or vitamin B12
        % Contains 2.25x more Dglc than RPMI
        % Also contains excess of basically everything
    elseif ismember({'DMEM'}, medium)
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEM');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end


    %% L15
        % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.80.html
        % Supports monkey kidney cells (HEP2) and primary explants of
        % embryonic and adult human tissue
        % includes galactose, phenol red, glutamine and sodium pyruvate
    elseif ismember({'L15'}, medium)
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'L15');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end
    %% McCoy
    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.83.html
    elseif ismember({'McCoy 5A'} , medium) % McCoy 5A
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'McCoy5A');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.76.html
    elseif ismember({'Iscove'} , medium) % IMDM
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'Iscove');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.126.html
    elseif ismember({'Waymouth'}, medium) % Waymouth
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'Waymouth');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.227.html
    elseif ismember({'DMEM:F12 (1:1)'}, medium) % 1:1 DMEM and F12
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEMF12');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.64.html
    elseif ismember({'HAM F12'}, medium) % F12
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'HAMF12');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/order/catalog/product/12571071?SID=srchsrp12571071
    % It is different from DMEM, but need to encode glutamine and nucleosides
    % into the code...
    elseif ismember({'alphaMEM'}, medium) % Alternative form to MEM
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'alphaMEM');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.114.html
    elseif ismember({'RPMI w Gln'}, medium) % RPMI ++ Lgln
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMIwGln');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.61.html
    elseif ismember({'HAM F10'}, medium) % F10
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'HAMF10');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % No one sells it  calc'd ratio
    elseif ismember({'DMEM:RPMI (2:1)'}, medium) % 2:1 DMEM:RPMI
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEMRPMI21');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % No one sells it  calc'd ratio
    elseif ismember({'MCDB105:M199'}, medium) % MCDB 105: Medium 199
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'MCDB105M199');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % https://www.thermofisher.com/us/en/home/technicalresources/mediaformulation.314.html
    elseif ismember({'Williams E Medium'}, medium) % Williams E Medium
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'Williams');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % NEED TO CHANGE IN FUTURE VERSION
    elseif ismember({'ACL4'}, medium) % ACL4
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEMF12');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % No one sells it  calc'd ratio
    elseif ismember({'RPMI:F12'}, medium) % 1:1 RPMI and F12
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMIF12');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % No one sells it  calc'd ratio
    elseif ismember({'DMEM:Iscove'}, medium) % 1:1 DMEM and Iscove 
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'DMEMIscove');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end

    % No one sells it  calc'd ratio  
    elseif ismember({'RPMI:Iscove'}, medium) % 1:1 RPMI and Iscove
        arr = readcell('./../../data/uptake.xlsx', 'Sheet', 'RPMIIscove');
        arr(1,:) = [];
        for rxn=1:length(arr)
            model.lb(find(ismember(model.rxns, arr(rxn, 2)))) = cell2mat(arr(rxn, 4));
        end
    end

end