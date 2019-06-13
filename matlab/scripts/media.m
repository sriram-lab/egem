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
    [num, txt] = xlsread('./../../data/uptake.xlsx', 'RPMI')
    for rxn=2:length()
    % Amino Acids
    model.lb(find(ismember(model.rxns, {'EX_gly(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_arg_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_asn_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_asp_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_cys_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_glu_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_gln_L(e)'}))) = -0.5;
    model.lb(find(ismember(model.rxns, {'EX_his_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_pro_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_leu_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_ile_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_lys_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_met_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_phe_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_ser_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_thr_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_trp_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_tyr_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_val_L(e)'}))) = -0.05;
    
    % Vitamins
    model.lb(find(ismember(model.rxns, {'EX_btn(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_chol(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_pnto_R(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_fol(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_ncam(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_pydx(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_aqcobal(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_inost(e)'}))) = -0.005;
    
    % Inorganic salts
    model.lb(find(ismember(model.rxns, {'EX_ca2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_k(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_na1(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cl(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pi(e)'}))) = -100;

    % Other
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;
    model.lb(find(ismember(model.rxns, {'EX_o2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_gthrd(e)'}))) = -0.05;
    
% If medium condition is not specified, set to RPMI
elseif ismember({'nan'}, medium)
    % Amino Acids
    model.lb(find(ismember(model.rxns, {'EX_gly(e)'}))) = -1;
    model.lb(find(ismember(model.rxns, {'EX_arg_L(e)'}))) = -0.10;
    model.lb(find(ismember(model.rxns, {'EX_asn_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_asp_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_cys_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_glu_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_gln_L(e)'}))) = -0.5;
    model.lb(find(ismember(model.rxns, {'EX_his_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_pro_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_leu_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_ile_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_lys_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_met_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_phe_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_ser_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_thr_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_trp_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_tyr_L(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_val_L(e)'}))) = -0.05;
    
    % Vitamins
    model.lb(find(ismember(model.rxns, {'EX_btn(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_chol(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_pnto_R(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_fol(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_ncam(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_pydx(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_aqcobal(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_inost(e)'}))) = -0.005;
    
    % Inorganic salts
    model.lb(find(ismember(model.rxns, {'EX_ca2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_k(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_na1(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cl(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pi(e)'}))) = -100;

    % Other
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;
    model.lb(find(ismember(model.rxns, {'EX_o2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_gthrd(e)'}))) = -0.05;

%% DMEM
    % https://www.thermofisher.com/order/catalog/product/12491015
    % No biotin or vitamin B12
    % Contains 2.25x more D-glc than RPMI
    % Also contains excess of basically everything
elseif ismember({'DMEM'}, medium)
    
    % Amino Acids
    model.lb(find(ismember(model.rxns, {'EX_gly(e)'}))) = -0.25; % 5x RPMI
    model.lb(find(ismember(model.rxns, {'EX_arg_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_asn_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_asp_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cys_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_glu_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_gln_L(e)'}))) = 0; % no glutamine
    model.lb(find(ismember(model.rxns, {'EX_his_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pro_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_leu_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_ile_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_lys_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_met_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_phe_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_ser_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_thr_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_trp_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_tyr_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_val_L(e)'}))) = -100;
    
    % Vitamins
    model.lb(find(ismember(model.rxns, {'EX_btn(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_chol(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pnto_R(e)'}))) = 0;
    model.lb(find(ismember(model.rxns, {'EX_fol(e)'}))) = -0.020;
    model.lb(find(ismember(model.rxns, {'EX_ncam(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pydx(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_aqcobal(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_inost(e)'}))) = -100;
    
    % Inorganic salts
    model.lb(find(ismember(model.rxns, {'EX_ca2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_k(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_na1(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cl(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pi(e)'}))) = -100;

    % Other
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*2.25;
    model.lb(find(ismember(model.rxns, {'EX_o2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_gthrd(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_pyr(e)'}))) = -100;

%% L15
    % https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.80.html
    % Supports monkey kidney cells (HEP-2) and primary explants of
    % embryonic and adult human tissue
    % includes galactose, phenol red, glutamine and sodium pyruvate
elseif ismember({'L15'}, medium) % No glc; some gal  
    
    % Amino Acids
    model.lb(find(ismember(model.rxns, {'EX_gly(e)'}))) = -0.25; % 5x RPMI
    model.lb(find(ismember(model.rxns, {'EX_ala_L(e)'}))) = -0.5;
    model.lb(find(ismember(model.rxns, {'EX_arg_L(e)'}))) = -0.10;
    model.lb(find(ismember(model.rxns, {'EX_asn_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_asp_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cys_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_glu_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_gln_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_his_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pro_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_leu_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_ile_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_lys_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_met_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_phe_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_ser_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_thr_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_trp_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_tyr_L(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_val_L(e)'}))) = -100;
    
    % Vitamins
    model.lb(find(ismember(model.rxns, {'EX_btn(e)'}))) = -0.005;
    model.lb(find(ismember(model.rxns, {'EX_chol(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pnto_R(e)'}))) = 0;
    model.lb(find(ismember(model.rxns, {'EX_fol(e)'}))) = -0.020;
    model.lb(find(ismember(model.rxns, {'EX_ncam(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pydx(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_aqcobal(e)'}))) = -0.05;
    model.lb(find(ismember(model.rxns, {'EX_inost(e)'}))) = -100;
    
    % Inorganic salts
    model.lb(find(ismember(model.rxns, {'EX_ca2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_k(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_na1(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_cl(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pi(e)'}))) = -100;

    % Other
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = 0; % No glc 
    model.lb(find(ismember(model.rxns, {'EX_gal(e)'}))) = -0.9; % Galactose
    model.lb(find(ismember(model.rxns, {'EX_o2(e)'}))) = -100;
    model.lb(find(ismember(model.rxns, {'EX_pyr(e)'}))) = -100; % excess pyruvate

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.83.html
elseif ismember({'McCoy 5A'} , medium) % McCoy 5A
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.76.html
elseif ismember({'Iscove'} , medium) % IMDM
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*4.5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.126.html
elseif ismember({'Waymouth'}, medium) % Waymouth
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.227.html
elseif ismember({'DMEM:F12 (1:1)'}, medium) % 1:1 DMEM and F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3.15/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.64.html
elseif ismember({'HAM F-12'}, medium) % F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.80/2;

% https://www.thermofisher.com/order/catalog/product/12571071?SID=srch-srp-12571071
% It is different from DMEM, but need to encode glutamine and nucleosides
% into the code...
elseif ismember({'alpha-MEM'}, medium) % Alternative form to MEM
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*4.5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
elseif ismember({'RPMI w Gln'}, medium) % RPMI ++ L-gln
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.61.html
elseif ismember({'HAM F-10'}, medium) % F10
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.1/2;

% No one sells it - calc'd ratio
elseif ismember({'DMEM:RPMI (2:1)'}, medium) % 2:1 DMEM:RPMI
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5/2;

% No one sells it - calc'd ratio
elseif ismember({'MCDB105:M199'}, medium) % MCDB 105: Medium 199
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -(5*0.72*1)/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.314.html
elseif ismember({'Williams E Medium'}, medium) % Williams E Medium
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% Same as DMEM:F12 
elseif ismember({'ACL-4'}, medium) % ACL-4
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3.15/2;  

% No one sells it - calc'd ratio
elseif ismember({'RPMI:F12'}, medium) % 1:1 RPMI and F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.8/2;   

% No one sells it - calc'd ratio
elseif ismember({'DMEM:Iscove'}, medium) % 1:1 DMEM and Iscove 
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% No one sells it - calc'd ratio  
elseif ismember({'RPMI:Iscove'}, medium) % 1:1 RPMI and Iscove
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*2.25;
end

end