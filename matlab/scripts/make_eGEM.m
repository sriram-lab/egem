%% Make eGEM
% @author: Scott Campit & Lauren Fane

% Initialize parameters
initCobraToolbox;
changeCobraSolver('gurobi');

% Load AcGEM eGEM (Shen et al., 2019)
load ./../shen-et-al/supplementary_software_code acetylation_model;
eGEM = acetylation_model;

% Check for duplicate reactions
eGEM = checkDuplicateRxn(eGEM,'S',1,1);

%% Reactions to remove:
eGEM = removeRxns(eGEM, {'peplyexn'});

%% Demand reactions that will be used in the study 
load './../vars/metabolites.mat'
<<<<<<< HEAD
compartment = char(metabolites(rxn, 4));
for m = 1:length(metabolites(:, 1))
    tmp_met = char(metabolites(m, 2));
    tmp = [tmp_met '[' compartment '] -> '];
    tmpname = char(metabolites(m, 1));
=======
compartment = 'n';
for m = 1:length(metabolites(:,1))
    tmp_met = char(metabolites(m,2));
    tmp = [tmp_met '[' compartment '] -> '];
    tmpname = char(metabolites(m,1));
>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967
    eGEM = addReaction(eGEM, tmpname, 'reactionFormula', tmp);
end

%% Histone position specific metabolic reactions
<<<<<<< HEAD
%[~, txt] = xlsread('metabolicmap.xlsx', 'Reactions');
=======
[~, txt] = xlsread('metabolicmap.xlsx', 'Reactions');

>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967

%% Add transport reactions to nucleus:
    
    % Methylation-related reactions
        % L-Methionine
        % ATP
    % Demethylation-related reactiosn
        % Succinate
        % Alpha-ketoglutarate
        % Fumarate
    % Acetylation-related reactions
        % Acetate
        % Pyruvate
        % Citrate (Shen et al., 2019)
    % One-carbon-related reactions
        % Formate
        % Glycine
        % Serine

eGEM = addReaction(eGEM, 'METtn',...
    'reactionFormula', 'met-L[c] <=> met-L[n]');
eGEM = addReaction(eGEM, 'SUCCtn',...
    'reactionFormula', 'succ[c] <=> succ[n]');
eGEM = addReaction(eGEM, 'AKGtn',...
    'reactionFormula', 'akg[c] <=> akg[n]');
eGEM = addReaction(eGEM, 'FUMtn',...
    'reactionFormula', 'fum[c] <=> fum[n]');
eGEM = addReaction(eGEM, 'FOLtn',...
    'reactionFormula', 'fol[c] <=> fol[n]');
eGEM = addReaction(eGEM, 'MALtn',...
    'reactionFormula', 'mal-L[c] <=> mal-L[n]');
%eGEM = addReaction(eGEM, 'GLYtn',...
%    'reactionFormula', 'gly[c] <=> fum[n]');
%eGEM = addReaction(eGEM, 'SERtn',...
%    'reactionFormula', 'ser-L[c] <=> ser-L[n]');
eGEM = addReaction(eGEM, 'ACtn',...
    'reactionFormula', 'ac[c] <=> ac[n]');
eGEM = addReaction(eGEM, 'PYRtn',...
    'reactionFormula', 'pyr[c] <=> pyr[n]');
%eGEM = addReaction(eGEM, 'FORtn',...
%    'reactionFormula', 'for[c] <=> for[n]');
eGEM = addReaction(eGEM, 'PEPtn',...
    'reactionFormula', 'pep[c] <=> pep[n]');

%% Add / modify existing Methionine cycle reactions

% Should try to add other methylation pathways, such as betaine, choline, etc
% as these pathways can also charge the methylation cycle. 

% MAT2A Production reaction
eGEM = addReaction(eGEM, 'METATn',...
    'reactionFormula', 'atp[n] + h2o[n] + met-L[n] -> amet[n] + pi[n] + ppi[n]',...
    'geneRule', '(MAT2A) and (MAT2B');
% HCy Production reaction
%eGEM = addReaction(eGEM, 'AHCn',...
%    'reactionFormula', 'ahcys[n] + h2o[n] -> adn[n] + hcys-L[n]',...
%    'geneRule', '(AHCY) or (AHCYL1) or (AHCYL1)');
% CBS
%eGEM = addReaction(eGEM, 'CYSTSn',...
%    'reactionFormula', 'hcys[n] + ser-L[n] -> cyst-L[n] + h2o[n]',...
%    'geneRule', '(CBS)');
% BHMT (somewhat - not doeGEMant nuclear mechanism)
%eGEM = addReaction(eGEM, 'BHMTn',...
%    'reactionFormula', 'glyb[n] + hcys__L[n] -> met-L[n] + dmgly[n]',...
%    'geneRule', '(BHMT) or (BHMT2)');

%% Improve specificity of histone methylation


%% Add demthylation reactions to nucleus (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732404/)
% Also added aldehyde detox pathway and source of Fe(II) and Fe(III)
    % Now that I have these as nodes in my eGEM, the challenge will be adding
    % some specificity in their performance

%LSD (Mono and di demethylation only)
eGEM = addReaction(eGEM, 'LSDDMT2n',...
   'reactionFormula', 'Ndmelys[n] + o2[n] + fad[n] -> fadh2[n] + h2o2[n] + mepepieGEMe[n]',...
   'geneRule', '(LSD1) or (LSD2)');
eGEM = addReaction(eGEM, 'PEPHYDROX2n',...
   'reactionFormula', 'mepepieGEMe[n] + h2o[n] -> Nmelys[n] + fald[n]',...
   'geneRule', '(LSD1) or (LSD2)');
eGEM = addReaction(eGEM, 'LSDDMT1n',...
   'reactionFormula', 'Nmelys[n]+ o2[n] + fad[n] -> fadh2[n] + h2o2[n] + pepieGEMe[n]',...
   'geneRule', '(LSD1) or (LSD2)');
eGEM = addReaction(eGEM, 'PEPHYDROX1n',...
   'reactionFormula', 'pepieGEMe[n] + h2o[n] -> peplys[n] + fald[n]',...
   'geneRule', '(LSD1) or (LSD2)');

%Dioxygenase (Mono-, di-, and tri-demethylation)
eGEM = addReaction(eGEM, 'JMJDMT3n',...
   'reactionFormula', 'Ntmelys[n] + akg[n] + o2[n] + fe2[n] -> Ndmelys[n] + fald[n] + co2[n] + fe3[n] + succ[n]',...
   'geneRule', '(JARID1A) or (JARID1B) or (JARID1C) or (JARID1D) or (JMJD3) or (JHDM3)');
eGEM = addReaction(eGEM, 'JMJDMT2n',...
   'reactionFormula', 'Ndmelys[n] + akg[n] + o2[n] + fe2[n] -> Nmelys[n] + fald[n] + co2[n] + fe3[n] + succ[n]',...
   'geneRule', 'JARID1A or JARD1B or JARID1C or JARID1D or JHDM2 or JHDM3 or PHF8 or JMJD3 or JHDM1');
eGEM = addReaction(eGEM, 'JMJDMT1n',...
   'reactionFormula', 'Nmelys[n] + akg[n] + o2[n] + fe2[n] -> peplys[n] + fald[n] + co2[n] + fe3[n] + succ[n]',...
   'geneRule', 'NO66 or JARID1B or JHDM2 or PHF8 or JHDM1');

% Aldehyde detoxification by ALDH to formate (GeneCards
%eGEM = addReaction(eGEM, 'ALDHn',...
%    'reactionFormula', 'fald[n] + nad[n] + h2o[n] -> for[n] + nadh[n] + h[n]',...
%    'geneRule', '(ALDH1A1) or (ALDH3A1) or (ALDH1B1) or (AKR1A1) or (ALDH1A3) or (ADH1B) or (ALDH7A1)');

% Iron reactions {Not sure how they are introduced so I created sink and demand reactions}
% % https://www.cell.com/trends/biochemical-sciences/pdf/S0968-0004(15)00237-6.pdf
<<<<<<< HEAD
eGEM = addReaction(eGEM, 'Fe2_sink',...
    'reactionFormula', 'fe2[n] <=> ');
eGEM = addReaction(eGEM, 'Fe3_sink',...
    'reactionFormula', 'fe3[n] <=> ');
eGEM = addReaction(eGEM, 'Fe2_demand',...
    'reactionFormula', 'fe2[n] -> ');
eGEM = addReaction(eGEM, 'Fe3_demand',...
    'reactionFormula', 'fe3[n] -> ');
eGEM = addReaction(eGEM, 'FE2tn', 'reactionFormula', 'fe2[c] -> fe3[n]');
=======
% eGEM = addReaction(eGEM, 'Fe2_sink',...
%     'reactionFormula', 'fe2[n] <=> ');
% eGEM = addReaction(eGEM, 'Fe3_sink',...
%     'reactionFormula', 'fe3[n] <=> ');
% eGEM = addReaction(eGEM, 'Fe2_demand',...
%     'reactionFormula', 'fe2[n] -> ');
% eGEM = addReaction(eGEM, 'Fe3_demand',...
%     'reactionFormula', 'fe3[n] -> ');
%eGEM = addReaction(eGEM, 'FE2tn', 'reactionFormula', 'fe2[c] -> fe3[n]');
>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967

%% FAD Pool (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4893201/)
% Since it is not known how flavin is introduced in the nucleus, I
% created a sink reaction and consumption of fad to fmn and amp.
    % RFK is only in cytoplasm (uniprot) -> no de novo synthesis in nucleus
% eGEM = addReaction(eGEM, 'FADn',...
%     'reactionFormula', 'fad[n] <=> '); 
% %eGEM = addReaction(eGEM, 'DM_FADn',...
% %    'reactionFormula', 'fad[n] -> '); 
% eGEM = addReaction(eGEM, 'FADDPn',...
%     'reactionFormula', 'fad[n] + h2o[n] -> amp[n] + fmn[n] + 2 h[n]',...
%     'geneRule', '(TKFC)');
% eGEM = addReaction(eGEM, 'FLAD1n',...
%     'reactionFormula', 'fmn[n] + atp[n] -> ppi[n] + fad[n]',...
%     'geneRule', '(FLAD1)');

%% Nuclear folate metabolism 
    
    % Nuclear folate metabolism is uncoupled to the methionine cycle, and
    % so adding these reactions may not be beneficial to eGEM acetylation
    % or methylation. 

    %(http://arjournals.annualreviews.org/doi/full/10.1146/annurev-nutr-071714-034441?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed)
    % https://academic.oup.com/advances/article/2/4/325/4591505
    % Sink reaction for folate: not known how folate is imported into
    % nucleus
    % Formate already has a transport reaction {FORtrn}
    % MTHFD1: Formate + THF + ATP + NADPH <=> 5,10meTHF + NADP+ + ADP + Pi
    % SHMT1, SHMT2a: THF + Ser <=> 5, 10meTHF + Gly
    % TYMS: THF+ dUMP -> DHF + dTMP
    % TK1 dT -> dTMP
        % Do I have to modify biomass reaction?
    % DHFR: THF + NADP+ <=> DHF + NADPH
%eGEM = addReaction(eGEM, 'Folate_sink',...
%    'reactionFormula', 'fol[n] <=> '); 
%eGEM = addReaction(eGEM, 'MTHFD1n',...
%    'reactionFormula', 'for[n] + thf[n] + atp[n] + nadph[n] <=> mlthf[n] + nadp[n] + adp[n] + pi[n]',...
%    'geneRule', '(MTHFD1)');
%eGEM = addReaction(eGEM, 'TYMS',...
%    'reactionFormula', 'mlthf[n] + dump[n] <=> dhf[n] + dtmp[n]',...
%    'geneRule', '(TYMS)');
%eGEM = addReaction(eGEM, 'TMDK1n',...
%    'reactionFormula', 'thymd[n] -> dtmp[n]',...
%    'geneRule', '(TMDK1)');
%eGEM = addReaction(eGEM, 'DHFRn',...
%    'reactionFormula', 'dhf[n] + nadph[n] -> thf[n] + nadp[n]',...
%    'geneRule', '(DHFR)');
%eGEM = addReaction(eGEM, 'SHMTn',...
%    'reactionFormula', 'thf[n] + ser[n] <=> gly[n] + mltf[n]',...
%    'geneRule', '(SHMT2)');

%% Central metabolism proteins

% While there is evidence that central metabolism enzymes co-localize to
% the nucleus for different functions other than metabolism, it is still
% undetereGEMed whether or not these metabolic enzymes actually perform
% metabolism, either due to a lack of complex parts that do no co-localize,
% PTMs, or other metabolic features. 

% https://www.nature.com/articles/s41580-018-0029-7
    % FH fum <=> mal
    % PKM2 pep + adp -> pyr + atp
    % PDC pyr <=> accoa
    % PFKFB4 f6p -> f26bp -> save for later
    % FBP1 f16bp + h2o -> f6p + pi
    % OGDH akg -> succ-coa + co2
    % GAPDH g3p + nad + p <=> 13bpg + nadh + h 
    % IMPDH imp + nad + h20 <=> xmp + nadh + h
    % GMPS atp + xmp + gln__L + h2o <=> amp + ppi + gmp + glu__L
    % ACSS2 ac + coa + atp -> accoa + adp + pi % the exact mechanism has not been detereGEMed,
    % but is an ATP-dependent reaction: https://www.cell.com/cancer-cell/pdf/S1535-6108(14)00511-X.pdf

% Not sure if these actually perform metabolic roles in nucleus - could be
% moonlighting
    
%eGEM = addReaction(eGEM, 'FUMn',...
%    'reactionFormula', 'fum[n] + h2o[n] <=> mal__L[n]',...
%    'geneRule', '(FH)');
%eGEM = addReaction(eGEM, 'FBPn',...
%    'reactionFormula', 'fdp[n] + h2o[n] -> f6p[n] + pi[n]',...
%    'geneRule', '(FBP1)');
%eGEM = addReaction(eGEM, 'AKGDn',...
%    'reactionFormula', 'akg[n] + coa[n] + nad[n] <=> co2[n] + nadh[n] + succoa[n]',...
%    'geneRule', '(((DLD and PDHX) and (OGDH and DLST)) or ((DLST and (DLD and PDHX) and (OGDH)))');
%eGEM = addReaction(eGEM, 'GAPDn',...
%    'reactionFormula', 'g3p[n] + nad[n] + pi[n] <=> 13dpg[n] + h[n] + nadh[n]',...
%    'geneRule', '(GAPDH)');
%eGEM = addReaction(eGEM, 'IMPDn',...
%    'reactionFormula', 'h2o[n] + imp[n] + nad[n] -> h[n] + nadh[n] + xmp[n]',...
%    'geneRule', '(IMPDH1) or (IMPDH2)');
%eGEM = addReaction(eGEM, 'GMPSn',...
%    'reactionFormula', 'atp[n] + nh3[n] + xmp[n] -> amp[n] + gmp[n] + 2.0 h[n] + ppi[n]',...
%    'geneRule', '(GMPS)');


%% Acetylation Reactions
eGEM = addReaction(eGEM, 'ACSS2n',...
    'reactionFormula', 'ac[n] + atp[n] + coa[n] <=> accoa[n] + amp[n] + ppi[n]',...
    'geneRule', '(ACSS2)');
eGEM = addReaction(eGEM, 'PDCn',...
    'reactionFormula', 'pyr[n] + coa[n] + nad[n] -> accoa[n] + co2[n] + nadh[n]',...
    'geneRule', '(DLAT and (DLD and PDHX) and (PDHA2 and PDHB)) or ((PDHA1 and PDHB) and (DLAT) and (DLD and PDHX))');
eGEM = addReaction(eGEM, 'PKM2n',...
    'reactionFormula', 'adp[n] + h[n] + pep[n] -> atp[n] + pyr[n]',...
    'geneRule', '(PKM2)');

%% Add more reactions - by LF
% %H3K9 methlation
% eGEM = addReaction(eGEM, 'Su(var)3-9',... 
%     'reactionFormula', '__lys[n] + SAM[n] <=> Nmelys[n] + SAH[n]', ... 
%     'checkDuplicate', 'true');
% %  Protein lysine + S-Adenosyl-L-methionine <=> Protein N6-methyl-L-lysine + S-Adenosyl-L-homocysteine
% %H3K4 methylation
% eGEM = addReaction(eGEM, 'H3K9MT',... 
%     'reactionFormula', '__l[n] + SAM[n] <=> Nmelys[n] + SAH[n]', ... 
%     'checkDuplicate', 'true');
%same chemical rxn as H3K9 methylation -> should i add it? How do i
%distinguish between locations?
%% save
eGEM = checkDuplicateRxn(eGEM,'S',1,1);
save('./../models/eGEM.mat', 'eGEM');