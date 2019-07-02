% Plot correlation values between AUC of growth inhibition of cell 
% lines and methylation flux

% Workflow: 1)Run make_eGEM  
cd '.\..\MeCorr';
% 2)Run methylVariables 3)Assign epsilon_methylation 4)Run last module of 
...MATLAB_CODE_methyl (3 sections)

% initCobraToolbox;
% changeCobraSolver('gurobi');

% s3 contains auc data. 
s3= readtable('Data\Ludlow2015_SmallMolecInformer.xlsx','Sheet','s3');
s3= table2array(s3);
% s2 matches cell line name to number assigned in this study
s2= readtable('Data\Ludlow2015_SmallMolecInformer.xlsx','Sheet','s2');

%The KMTi's and KDMi's in Ludlow data s3
%index_KMTi= [3, 178, 280, 374, 380, 421, 431, 475];
%hmei_list= {'BRD-A02303741','BIX-01294','methylstat','QW-BI-011',...
%    'UNC0321','CBB-1007','UNC0638','GSK-J4'};
% drugs 280, 421, and 475 are KDMi's. Others are KMTi's
%index_DNMTi = [2 , 78, 181, 402, 444];     % DNMT = DNA MT

% Get all 3 rows of s3 for drug_auc_me. Replaces drug_auc_expt
% 7/1: drug_auc_expt has all drugs tested. drug_auc_me only has methyl drugs 
% advantage: less to search. matlab_code already extracts interested drugs
% (variable hmei_auc_expt_me).
drug_auc_me= s3;

% drug_auc_me= zeros(664*8, 3);
% drug_auc_me(1:597, :)= s3(1293:1889, :); % drug 3
% drug_auc_me(665:1291, :)= s3(96133:96759, :); %drug 178
% drug_auc_me(1328:1766, :)= s3(149039:149477, :); %drug 280
% drug_auc_me(1992:2435, :)= s3(201699:202142, :); %drug 374
% drug_auc_me(2656:3302, :)= s3(204605:205251, :); %drug 380
% drug_auc_me(3320:3630, :)= s3(228209:228519, :); %drug 421
% drug_auc_me(3984:4606, :)= s3(232799:233421, :); %drug 431
% drug_auc_me(4648:4718, :)= s3(257167:257237, :); %drug 475

% Col 1 is drug id #, col 2 is cell line id #, col 3 is auc. Switch col 2 
...and 3 to match drug_auc_expt
celllineid_tmp= drug_auc_me(:, 2);
drug_auc_me(:, 2)= drug_auc_me(:, 3);
drug_auc_me(:, 3)= celllineid_tmp;
drug_auc_me= array2table(drug_auc_me);
% Rename drug_auc_me Var Names  so they match that of
% ctd2compoundidname_id_me and describe the columns
drug_auc_me.Properties.VariableNames{'drug_auc_me1'}='index_cpd';
drug_auc_me.Properties.VariableNames{'drug_auc_me2'}='auc';
drug_auc_me.Properties.VariableNames{'drug_auc_me3'}='index_ccl';

% Other variables to replace with methylation data
ctd2clidname_id_me= s2(:, 1);

ctd2clidname_me= [s2(:,2), s2(:,4), s2(:,5)]; 
ctd2clidname_me.Properties.VariableNames{'cell_line_name'}='celllinenames_ccle1';

s1= readtable('Data\Ludlow2015_SmallMolecInformer.xlsx','Sheet','s1','Range','A:B');
ctd2cpdidname_id_me= s1(:, 1); 
ctd2cpdidname_name_me= s1(:, 2);

% Below: written into MATALB_CODE_methyl
% exptidcelllinemediamatch= array2table(exptidcelllinemediamatch);
% celllinenames_ccle1= cell2table(celllinenames_ccle1);

% Replace hcommon_exptdat with KMTi_auc
KMTi_auc= zeros(664, 8); % Ludlow tested 664 cell lines. 8 methyl drugs

KMTi_auc(1:597, 1)= s3(1293:1889, 3); % drug 3
KMTi_auc(1:627, 2)= s3(96133:96759, 3); %drug 178
KMTi_auc(1:439, 3)= s3(149039:149477, 3); %drug 280
KMTi_auc(1:444, 4)= s3(201699:202142, 3); %drug 374
KMTi_auc(1:647, 5)= s3(204605:205251, 3); %drug 380
KMTi_auc(1:311, 6)= s3(228209:228519, 3); %drug 421
KMTi_auc(1:623, 7)= s3(232799:233421, 3); %drug 431
KMTi_auc(1:71, 8)= s3(257167:257237, 3); %drug 475