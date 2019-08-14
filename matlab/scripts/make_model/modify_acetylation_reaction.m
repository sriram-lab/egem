% modify_acetylation_reaction

model = readCbModel('RECON1.xml');
load ./../shen-et-al/supplementary_software_code acetylation_model;

acetylation_model = removeRxns(acetylation_model, 'EX_retpalm');

% Need to modify these variables to also update other fields
acetylation_model.metCharge = model.metCharge;
acetylation_model.csense = model.csense;
acetylation_model.metNames = model.metNames;
acetylation_model.metHMDBID = model.metHMDBID;
acetylation_model.metKEGGID = model.metKEGGID;
acetylation_model.metChEBIID = model.metChEBIID;
acetylation_model.metMetaNetXID = model.metMetaNetXID;
acetylation_model.rxnECNumbers = model.rxnECNumbers;
acetylation_model.rxnKEGGID = model.rxnKEGGID;
acetylation_model.rxnMetaNetXID = model.rxnMetaNetXID;
acetylation_model.rxnSBOTerms = model.rxnSBOTerms;

clear all;