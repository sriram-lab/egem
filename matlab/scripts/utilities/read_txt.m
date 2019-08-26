%% @author: Scott Campit 
%function [num, txt] = read_txt(filePath)
%    %var = textread(filePath, '%s%s%s', 'whitespace', '\n');
%    [num, txt] = xlsread(filePath);
%end

<<<<<<< HEAD
eGEM.mets = strrep(eGEM.mets, '[', '_');
eGEM.mets = strrep(eGEM.mets, ']', '');

% Will not read because there is error in stoichiometric matrix
biomassRxn = {'0.505626 ala__L_c + 0.35926 arg__L_c + 0.279425 asn__L_c + 0.352607 asp__L_c + 20.704451 atp_c + 0.039036 ctp_c + 0.046571 cys__L_c + 0.275194 g6p_c + 0.325996 gln__L_c + 0.385872 glu__L_c + 0.538891 gly_c + 0.036117 gtp_c + 20.650823 h2o_c + 0.126406 his__L_c + 0.286078 ile__L_c + 0.545544 leu__L_c + 0.592114 lys__L_c + 0.153018 met__L_c + 0.259466 phe__L_c + 0.412484 pro__L_c + 0.392525 ser__L_c + 0.31269 thr__L_c + 0.013306 trp__L_c + 0.159671 tyr__L_c + 0.053446 utp_c + 0.352607 val__L_c + 0.020401 chsterol_c + 0.011658 clpn_hs_c + 0.023315 pail_hs_c + 0.154463 pchol_hs_c + 0.055374 pe_hs_c + 0.002914 pglyc_hs_c + 0.005829 ps_hs_c + 0.017486 sphmyln_hs_c + 0.013183 datp_n + 0.009442 dctp_n + 0.009898 dgtp_n + 0.013091 dttp_n ? 20.650823 adp_c + 20.650823 h_c + 20.650823 pi_c'};
new_model = addReaction(recon1, 'biomassReaction', ...
    'reactionFormula', biomassRxn, 'rev', 0, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Biomass and maintenance functions');

% This reaction has the biomass reaction I want and all the reactions.
model = readCbModel('Recon3D.xml');
=======
>>>>>>> 931f91b80f9523ac8900fe7f7422d0229f08fa82
