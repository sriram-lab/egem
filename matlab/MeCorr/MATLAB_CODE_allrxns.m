      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % this code reproduces all the main text figures from the manuscript (Shen et al)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk methylation model:
% 1) cd ./../scripts. Run make_eGEM to create bulk methylation model
% 2) cd ./../MeCorr. Run section 1 and 3 of this script. Only run Section 2 if you want to 
% run the long for-loop
% 3) Parameters of interest to vary: epsilon_methylation, model2
model2 = min_model;
epsilon_methylation = 1E-1;
rxnpos  = find(ismember(min_model.rxns,'LYSMTF1n'));
minfluxflag = 0; % 0: no Pfba, 1: Pfba 
% 4) Change name of files to which variables are saved. End of section 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % impact of basal metabolic state of CCLE cell lines on sensitivity to demethylase inhibitors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load supplementary_software_code  ctd2celllineidname exptidcelllinemediamatch ctd2celllineidname_id r* 
%data from seashore-ludlow study, contains cell line names, growth media
load supplementary_software_code  ctd2compoundidname_id drug_auc_expt ctd2compoundidname_name 
%data from seashore-ludlow study, contains drug names, drug sensitivity data
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz  
%contains CCLE cell line names, gene expression data (z-transformed)
load supplementary_software_code hcommon_exptdat hcommon1 
%data from seashore-ludlow study, contains drug names, drug sensitivity data for cell lines that were screened against all 8 hme inhibitors
load supplementary_software_code hdacexpfcs hdacexpallgeneids 
%contains gene expression data after treatment with hdac inhibitors

% Converting arrays to tables and renaming some Variables to make
% caompatible with function ismember
ctd2celllineidname= array2table(ctd2celllineidname);
ctd2celllineidname.Properties.VariableNames{'ctd2celllineidname1'}= 'ccl_name';
ctd2celllineidname_id= array2table(ctd2celllineidname_id);
ctd2celllineidname_id.Properties.VariableNames{'ctd2celllineidname_id'}= 'master_ccl_id';
ctd2compoundidname_name= array2table(ctd2compoundidname_name);
ctd2compoundidname_name.Properties.VariableNames{'ctd2compoundidname_name'}= 'cpd_name';
ctd2compoundidname_id= array2table(ctd2compoundidname_id);
ctd2compoundidname_id.Properties.VariableNames{'ctd2compoundidname_id'}= 'cpd_id';
celllinenames_ccle1= array2table(celllinenames_ccle1);
celllinenames_ccle1.Properties.VariableNames{'celllinenames_ccle1'}= 'ccl_name';
drug_auc_expt= array2table(drug_auc_expt);
drug_auc_expt.Properties.VariableNames{'drug_auc_expt3'}= 'cpd_id';
drug_auc_expt.Properties.VariableNames{'drug_auc_expt1'}= 'experiment_id';
drug_auc_expt.Properties.VariableNames{'drug_auc_expt2'}= 'auc';
exptidcelllinemediamatch=exptidcelllinemediamatch(2:end, :);
exptidcelllinemediamatch= array2table(exptidcelllinemediamatch); %s4
exptidcelllinemediamatch.Properties.VariableNames{'exptidcelllinemediamatch1'}= 'experiment_id';
exptidcelllinemediamatch.Properties.VariableNames{'exptidcelllinemediamatch2'}= 'master_ccl_id';
hcommon1= array2table(hcommon1);
hcommon1.Properties.VariableNames{'hcommon1'}= 'experiment_id';

MODE = 1;  % reaction (1) or gene list (0)
epsilon = 1E-2; rho = 1;
kappa = 1; % parameters for integrating transcriptomics data. kappa is the strength of down regulation of genes (Chandrasekaran & Price, PNAS, 2010)
hmetransint = 1;
basalflag = 1;

hmei_list = {'BRD-A02303741';'BIX-01294';'methylstat';'QW-BI-011';...
    'UNC0321';'CBB-1007';'UNC0638';'GSK-J4'};
hmei_list= cell2table(hmei_list);
hmei_list.Properties.VariableNames{'hmei_list'}='cpd_name';

% If I don't want to run next for-loop:
%load('VariablesSaved\fluxstate_gurobi');
%load('VariablesSaved\grate_ccle_exp_soft');
%% long for-loop (1:1035)
grate_ccle_exp_soft = NaN(height(exptidcelllinemediamatch),3); %contains acetylation flux based on basal metabolic state
fluxes_allrxns= NaN(height(exptidcelllinemediamatch),length(model2.rxns));
g_rate= NaN(height(exptidcelllinemediamatch), length(model2.rxns));
for i = 1:height(exptidcelllinemediamatch)
    % match cell line data in CCLE with CTD2
    ii = find(ismember(ctd2celllineidname_id, exptidcelllinemediamatch(i,2)));
    iii = find(ismember(celllinenames_ccle1, ctd2celllineidname(ii,1)));
    if ~isempty(iii)
        iii  = iii(1);
%        model2 = min_model; 
%        model2 = metabolicmodel;
      
         %find up and down-regulated genes in each cell line
        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) > 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) < -2));
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  set the glucose uptake based on media
        % default glucose is -5 for rpmi
        if r1(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5;% no change rpmi
        elseif r2(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*1/2;% dmem with low glucose
        elseif r3(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*1/2;% Emem..
        elseif r4(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*3/2;% mccoy
        elseif r5(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*1/2;% mem
        elseif r6(i)
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*3.15/2;% dmem:f12
        end
        
        %find reactions from differentially expressed genes
        onreactions= findRxnsFromGenes(model2, ongenes);
        onreactions= struct2cell(onreactions);
        for jj=1:length(onreactions)
            onreactions(jj)= onreactions{jj}(1);
        end
        offreactions= findRxnsFromGenes(model2, offgenes);
        offreactions= struct2cell(offreactions);
        for jj=1:length(offreactions)
            offreactions(jj)= offreactions{jj}(1);
        end
        % Below 2 lines work for acetyl model, not min methyl model.
%         [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
%         [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);        
        
        disp(i)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get basal metabolic state based on transcriptome
%         if basalflag
%             [fluxes, grate, solverobj_ccle] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[], minfluxflag); % impact on growth
%             grate_ccle_exp_soft(i,1:2)= grate; % first 2 columns are basal met flux. 
%             ...3rd column will be met flux with a rxn maximized
%             fluxes_allrxns(i,:) = fluxes; % correlate each row of fluxes_allreactions with auc
%             model2.c(rxnpos) = epsilon_methylation ;
%             [fluxstate_gurobi] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[],minfluxflag);
%             grate_ccle_exp_soft(i,3) = fluxstate_gurobi(rxnpos); %methylation flux
%         end
        
        for rxncount = 1:length(model2.rxns)
            if basalflag
                %minfluxflag = 1;
                %[~, grate, solverobj_ccle] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[], minfluxflag); % impact on growth
                model2.c(rxncount) = epsilon_methylation ;
                [fluxstate_gurobi] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[],minfluxflag);
                g_rate(i,rxncount) = fluxstate_gurobi(rxncount); %methylation flux
            end
            disp(rxncount)
        end
    else
        grate_ccle_exp_soft(i,:) = NaN;
    end
end

% figure; h = histogram(grate_ccle_exp_soft(:,3),70);
%  xlabel('Predicted methylation flux')
%  ylabel('Total cell lines')
% title('Distribution of methylation flux among CCLE cell lines','fontweight','normal')

% save('VariablesSaved\g_rate_1E-1', 'g_rate');
% save('VariablesSaved\fluxesAll_1E-6', 'fluxes_allrxns');
% save('VariablesSaved\grate_1E-6', 'grate_ccle_exp_soft');
% save('VariablesSaved\fluxstate_1E-6', 'fluxstate_gurobi');
%% Correlation between flux and auc for a reaction
% flux_allreactions: 1031 cell lines by 3777 reactions
% drug_auc_expt: use auc values for the 8 methyl drugs
% 8 x 3000 rhos
rho= NaN(height(hmei_list),length(model2.rxns));
rho2= NaN(height(hmei_list),length(model2.rxns));

for j = 1:height(hmei_list)
    fx = find(ismember(ctd2compoundidname_name, hmei_list(j,1)));
    ix = ismember(drug_auc_expt(:,3), ctd2compoundidname_id(fx,1)); sum(ix);
    hmei_auc_dat = drug_auc_expt(ix,:);
    
    % Only uses the experiments of fluxes_allrxns that also are in hmei_auc_dat
    [ix, pos] = ismember(hmei_auc_dat(:,1), exptidcelllinemediamatch(:,1)); sum(ix); 
    v1 = hmei_auc_dat(:,2); v1= table2array(v1);
    for nCol=1:size(fluxes_allrxns,2)
        v4= fluxes_allrxns(pos,nCol);
        rho(j,nCol)= corr(v4, v1, 'rows', 'complete');
    end
%     for nCol=1:size(fluxes_allrxns,2)
%         v4= fluxes_allrxns(pos,nCol);
%         ix2= find(~isnan(v4));
%         v1= v1(ix2);
%         v4= v4(~isnan(v4));
%         rho2(j,nCol)= corr2(v4, v1);
%     end
end
threshold= [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2];
j2=7;
sigIndTF= (abs(rho) > threshold(j2));
nSigExpt= sum(sigIndTF); % sum number of signif expts per rxn (sum each column)
sigRxnInd= (nSigExpt >= 1);
sum(sigRxnInd)
sigRxn= model2.rxns(sigRxnInd);
disp(sigRxn) 
%% Calculate correlation between flux and auc. Each reaction maximized.
rho3= NaN(height(hmei_list),length(model2.rxns));
for j= 1:height(hmei_list)
    % Extract the auc for a drug in drug_auc_expt
    fx = find(ismember(ctd2compoundidname_name, hmei_list(j,1)));
    ix = ismember(drug_auc_expt(:,3), ctd2compoundidname_id(fx,1));
    hmei_auc_dat = drug_auc_expt(ix,:);
    % Only uses the experiments of fluxes_allrxns that also are in hmei_auc_dat
    [ix, pos] = ismember(hmei_auc_dat(:,1), exptidcelllinemediamatch(:,1));
    v1 = hmei_auc_dat(:,2); 
    v1= table2array(v1);
    
    for nCol= 1:size(g_rate,2)
        v5= g_rate(pos, nCol);
        rho3(j,nCol)= corr(v5, v1, 'rows', 'complete');
    end
end
save('VariablesSaved\rho3_1E-1', rho3);

sigRho3Ind= false(1,length(model2.rxns));
for j2= 13:13%length(model2.rxns)
    ix= ~isnan(rho3(:,j2));
    sigRho3Ind(1,j2)= (sum(ix)>0);
end
sigRxn3= model2.rxns(sigRho3Ind);
disp(sigRxn3)
save('VariablesSabed\sigRxn3_1E-1', sigRxn3);
%% predicting sensitivity to hme inhibitors based on basal metabolic state - comparison with drug sensitivity data from seashore-ludlow study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars r_basal
r_basal(height(hmei_list))= struct('rho',0, 'p',0, 'rlowerbound',0, 'rupperbound',0);
pp_flux= zeros(height(hmei_list), 3);
pp_grate= zeros(height(hmei_list), 3);
for j = 1:height(hmei_list)
    fx = find(ismember(ctd2compoundidname_name, hmei_list(j,1)));
    ix = ismember(drug_auc_expt(:,3), ctd2compoundidname_id(fx,1));
    sum(ix) % 847. ix is a 1D array. Each drug appears multiple times (appears for each medium in drug_auc_expt)
    hme_auc_dat = drug_auc_expt(ix,:);
    
    [ix, pos] = ismember(hme_auc_dat(:,1), exptidcelllinemediamatch(:,1)); sum(ix);
    v1 = hme_auc_dat(:,2); v1= table2array(v1);
    v2 = grate_ccle_exp_soft(pos,:); % growth rate
    v3 = fluxstate_gurobi(pos,:);
    
%     % Scatter plots of methylation (metabolic) flux & calculate correlation
%     figure;
%     scatter(v1, v3)
%     ylabel({'Methylation Flux'},'fontname','helvetica'); %'fontweight','bold')
%     xlabel({'Sensitivity (AUC)'; hmei_list.cpd_name{j}},'fontname','helvetica');
%     ylim([0, 0.06])
%     [R, P, RL, RU]= corrcoef(v1, v3);
%     r_basal(j).rho=R; r_basal(j).p=P; r_basal(j).rlowerbound=RL; r_basal(j).rupperbound=RU;
%     str=sprintf('r= %1.3f',r_basal(j).rho(1,2));
%     T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%     set(T, 'fontsize', 14, 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
    
    % Boxplot pairs of methylation flux (Actually graphing metabolic flux, 
    ...but we hypothesize met flux approximates methyl flux)
    ix0 = (v3 < 0.05);sum(ix0); % 0.05 is  the default threshold value
    ix01 = (v3 > 0.05);sum(ix01);
    [hh, pp_flux(j,3)] = ttest2(v1(ix0), v1(ix01)); %
    
    groups = NaN(size(ix0));
    groups(ix0) = 1; groups(ix01) = 2;
    vv = NaN(length(v1), 2);
    vv(1:sum(groups == 1),1) = v1(groups == 1);
    vv(1:sum(groups == 2),2) = v1(groups == 2);
%     
%     figure;
%     %clf; UnivarScatter(vv,'Width', 0.3, 'PointSize', 11,'MarkerEdgeColor','w','LineWidth',0.1);%, 'markerfacealpha',0.5);
%     %hold on;
%     bh = boxplot(v1, groups,'symbol','');
%     set(gca,'xticklabel',{'Low flux','High flux'})
%     ylabel({'Sensitivity (AUC)'},'fontname','helvetica');%,'fontweight','bold');
%     xlabel(hmei_list{j,1});
%     %ylim([0 20])
%     str=sprintf('p= %1.4f',pp_flux(j,3));
%     T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%     set(T, 'fontsize', 14, 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
    
    ix0 = (v3 <= prctile(v3, 25)); sum(ix0);
    ix01 = (v3 > prctile(v3, 25)); sum(ix01);
    [hh, pp_flux(j,1)] = ttest2(v1(ix0), v1(ix01));
    
    ix0 = (v3 <= prctile(v3, 50)); sum(ix0);
    ix01 = (v3 > prctile(v3, 50)); sum(ix01);
    [hh, pp_flux(j,2)] = ttest2(v1(ix0), v1(ix01));
    
    % Boxplot pairs of growth rate %%%%%%
    ix0 = (v2(:,2) < 0.05);sum(ix0); % 0.05 is  the default threshold value
    ix01 = (v2(:,2) > 0.05);sum(ix01);
    [hh, pp_grate(j,3)] = ttest2(v1(ix0), v1(ix01)); %
    
    groups = NaN(size(ix0));
    groups(ix0) = 1; groups(ix01) = 2;
    vv = NaN(length(v1), 2);
    vv(1:sum(groups == 1),1) = v1(groups == 1);
    vv(1:sum(groups == 2),2) = v1(groups == 2);
    
%     figure;
%     %clf; UnivarScatter(vv,'Width', 0.3, 'PointSize', 11,'MarkerEdgeColor','w','LineWidth',0.1);%, 'markerfacealpha',0.5);
%     %hold on;
%     bh = boxplot(v1, groups,'symbol','');
%     set(gca,'xticklabel',{'Low growth','High growth'})
%     ylabel({'Sensitivity (AUC)'},'fontname','helvetica');%,'fontweight','bold');
%     xlabel(hmei_list{j,1});
%     %ylim([0 20])
%     str=sprintf('p= %1.4f',pp_grate(j,3));
%     T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%     set(T, 'fontsize', 14, 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
    
    ix0 = (v2(:,3) <= prctile(v2(:,3), 25)); sum(ix0);
    ix01 = (v2(:,3) > prctile(v2(:,3), 25)); sum(ix01);
    [hh, pp_grate(j,1)] = ttest2(v1(ix0), v1(ix01));
    
    ix0 = (v2(:,3) <= prctile(v2(:,3), 50)); sum(ix0);
    ix01 = (v2(:,3) > prctile(v2(:,3), 50)); sum(ix01);
    [hh, pp_grate(j,2)] = ttest2(v1(ix0), v1(ix01));

%         ix0 = (v2(:,2) <= prctile(v2(:,2), 25)); sum(ix0);
%         ix01 = (v2(:,2) > prctile(v2(:,2), 25)); sum(ix01);
%         [hh, pp_grate(j,1)] = ttest2(v1(ix0), v1(ix01));
% 
%         ix0 = (v2(:,2) <= prctile(v2(:,2), 50)); sum(ix0);
%         ix01 = (v2(:,2) > prctile(v2(:,2), 50)); sum(ix01);
%         [hh, pp_grate(j,2)] = ttest2(v1(ix0), v1(ix01));
end
r_basal= struct2table(r_basal); 
r_basal.Properties.RowNames= {hmei_list.cpd_name{1:8}};
disp(r_basal)
pp_flux= array2table(pp_flux);
pp_flux.Properties.RowNames= {hmei_list.cpd_name{1:8}};
pp_flux.Properties.VariableNames{'pp_flux3'}= 'pf_hilo';
pp_flux.Properties.VariableNames{'pp_flux1'}= 'pf_25prctile';
pp_flux.Properties.VariableNames{'pp_flux2'}= 'pf_50prctile';
disp(pp_flux) %t-test p-values for each drug
%save('VariablesSaved\pp_flux_1E-2_pFBA', 'pp_flux');
%load('VariablesSaved\pp_flux_1E-2_pFBA');
pp_grate= array2table(pp_grate);
pp_grate.Properties.RowNames= {hmei_list.cpd_name{1:8}};
pp_grate.Properties.VariableNames{'pp_grate3'}= 'pg_hilo';
pp_grate.Properties.VariableNames{'pp_grate1'}= 'pg_25prctile';
pp_grate.Properties.VariableNames{'pp_grate2'}= 'pg_50prctile';
disp(pp_grate)
%save('VariablesSaved\pp_grate_1E-6_pFBA', 'pp_grate');