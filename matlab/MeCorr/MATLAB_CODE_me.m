      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % this code reproduces all the main text figures from the manuscript (Shen et al)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk methylation model:
% 1) cd ./../scripts. Run make_eGEM to create bulk methylation model
model2 = min_model;
% 2) cd ./../MeCorr. Run section 1 and 3 of this script. Only run Section 2 if you want to 
% run the long for-loop
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
exptidcelllinemediamatch= array2table(exptidcelllinemediamatch); %s4
exptidcelllinemediamatch.Properties.VariableNames{'exptidcelllinemediamatch1'}= 'experiment_id';
exptidcelllinemediamatch.Properties.VariableNames{'exptidcelllinemediamatch2'}= 'master_ccl_id';
hcommon1= array2table(hcommon1);
hcommon1.Properties.VariableNames{'hcommon1'}= 'experiment_id';

rxnpos  = find(ismember(min_model.rxns,'LYSMTF1n')); % rxnpos = 2451
epsilon_methylation = 1E-2; % or 1E-1

MODE = 1;  % reaction (1) or gene list (0)
epsilon = 1E-2; rho = 1;
kappa = 1; % parameters for integrating transcriptomics data. kappa is the strength of down regulation of genes (Chandrasekaran & Price, PNAS, 2010)
minfluxflag = 0; % 0: no Pfba, 1: Pfba  
hmetransint = 1;
basalflag = 1;

% If I don't want to run next for-loop:
%load('VariablesSaved\fluxstate_gurobi');
%load('VariablesSaved\grate_ccle_exp_soft');
%% long for-loop (1:1035)
grate_ccle_exp_soft = NaN(height(exptidcelllinemediamatch),2); %contains acetylation flux based on basal metabolic state
%grate_ccle_exp_soft_hdacsign = NaN(height(exptidcelllinemediamatch),8); %contains acetylation flux based on basal metabolic state and impact of each hme inhibitor teatment

for i = 1:height(exptidcelllinemediamatch)
    % match cell line data in CCLE with CTD2
    ii = find(ismember(ctd2celllineidname_id, exptidcelllinemediamatch(i,2)));
    iii = find(ismember(celllinenames_ccle1, ctd2celllineidname(ii,1)));
    if ~isempty(iii)
        iii  = iii(1);
%         model2 = min_model; 
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
        if basalflag
            [~, grate_ccle_exp_soft, solverobj_ccle] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[], minfluxflag); % impact on growth
            model2.c(rxnpos) = epsilon_methylation ;
            [fluxstate_gurobi] = constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[],minfluxflag);
            grate_ccle_exp_soft(i,2) = fluxstate_gurobi(rxnpos); %acetylation flux
        end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % hme inhibitor impact on transcriptome
%         if hmetransint
%             for kk = 4:-1:1   % match drug order
%                 ongenesh = hdacexpallgeneids((hdacexpfcs(:,kk) > 1.1) & (hdacexpfcs(:,kk + 4) < 0.01));
%                 ongenesh = intersect(ongenesh, model2.genes);
%                 offgenesh = hdacexpallgeneids((hdacexpfcs(:,kk) < 0.9) & (hdacexpfcs(:,kk + 4) < 0.01));
%                 offgenesh = intersect(offgenesh, model2.genes);
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 [~,~,onreactionsh,~] =  deleteModelGenes(model2, ongenesh);
%                 [~,~,offreactionsh,~] =  deleteModelGenes(model2, offgenesh);
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % remove conflicting transcripts
%                 onreactionsh0 = setdiff(onreactionsh, offreactions);
%                 offreactionsh0 = setdiff(offreactionsh, onreactions);
%                 onreactions0 = setdiff(onreactions, offreactionsh);
%                 offreactions0 = setdiff(offreactions, onreactionsh);
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 onreactions1h = [onreactions0;  onreactionsh0];
%                 offreactions1h = [offreactions0;  offreactionsh0];
%                 epsilon2 = [repmat(0, size(offreactions1h))];
%                 model2.c(rxnpos) = 0;
%          %        [fluxstate_gurobi,grate_ccle_exp_soft_hdacsign(i,kk)] =  constrain_flux_regulation(model2,onreactions1h,offreactions1h,kappa,rho,epsilon,MODE, epsilon2,minfluxflag); % impact on growth
%                 model2.c(rxnpos) = epsilon_methylation ;
%                   [fluxstate_gurobi] =  constrain_flux_regulation(model2,onreactions1h,offreactions1h,kappa,rho,epsilon,MODE,epsilon2,minfluxflag);
%                  grate_ccle_exp_soft_hdacsign(i,kk + 4) = fluxstate_gurobi(rxnpos); %acetylation flux
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             end
%         end

    else
        grate_ccle_exp_soft(i,:) = NaN;
%        grate_ccle_exp_soft_hdacsign(i,:) = NaN;
    end
end


figure; h = histogram(grate_ccle_exp_soft(:,2),70);
 xlabel('Predicted methylation flux')
 ylabel('Total cell lines')
title('Distribution of methylation flux among CCLE cell lines','fontweight','normal')

save('VariablesSaved\grate_ccle_exp_soft', 'grate_ccle_exp_soft');
save('VariablesSaved\fluxstate_gurobi', 'fluxstate_gurobi');
save('VariablesSaved\solverobj_ccle', 'solverobj_ccle');
%% predicting sensitivity to hme inhibitors based on basal metabolic state - comparison with drug sensitivity data from seashore-ludlow study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hmei_list = {'BRD-A02303741';'BIX-01294';'methylstat';'QW-BI-011';...
    'UNC0321';'CBB-1007';'UNC0638';'GSK-J4'};
hmei_list= cell2table(hmei_list);
hmei_list.Properties.VariableNames{'hmei_list'}='cpd_name';
r_basal= cell(height(hmei_list),1);
%r_basal(height(hmei_list))= struct('r',0, 'p',0, 'RL',0, 'RP',0);
pp_basal= zeros(height(hmei_list), 3);
for j = 1:height(hmei_list)
    fx = find(ismember(ctd2compoundidname_name, hmei_list(j,1)));
    ix = ismember(drug_auc_expt(:,3), ctd2compoundidname_id(fx,1));
    sum(ix) % 847. ix is a 1D array. Each drug appears multiple times (appears for each medium in drug_auc_expt)
    hme_auc_dat = drug_auc_expt(ix,:);
    
    [ix, pos] = ismember(hme_auc_dat(:,1), exptidcelllinemediamatch(:,1)); sum(ix); %
    v1 = hme_auc_dat(:,2); v1= table2array(v1);
    %v2 = grate_ccle_exp_soft(pos,:);
    v3 = fluxstate_gurobi(pos,:);
    
    figure;
    scatter(v1, v3)
    ylabel({'Methylation Flux'},'fontname','helvetica'); %'fontweight','bold')
    xlabel({'Sensitivity (AUC)'; hmei_list.cpd_name{j}},'fontname','helvetica');
    ylim([0, 0.06])
    %r_basal(j).r=R; r_basal(j).p=P; r_basal(j).rl=RL; r_basal(j).ru=RU;
    r_basal{j}= corrcoef(v1, v3);
    str=sprintf('r= %1.3f',r_basal{j}(1,2));
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
%     ix0 = (v2(:,2) < 0.05);sum(ix0); % 0.05 is  the default threshold value
%     ix01 = (v2(:,2) > 0.05);sum(ix01);
%     [hh, pp_basal(j,3)] = ttest2(v1(ix0), v1(ix01)); %
%     
%     groups = NaN(size(ix0));
%     groups(ix0) = 1; groups(ix01) = 2;
%     vv = NaN(length(v1), 2);
%     vv(1:sum(groups == 1),1) = v1(groups == 1);
%     vv(1:sum(groups == 2),2) = v1(groups == 2);
%     
%     figure;
%     %clf; UnivarScatter(vv,'Width', 0.3, 'PointSize', 11,'MarkerEdgeColor','w','LineWidth',0.1);%, 'markerfacealpha',0.5);
%     %hold on;
%     bh = boxplot(v1, groups,'symbol','')
%     set(gca,'xticklabel',{'Low flux','High flux'})
%     ylabel({'Sensitivity (AUC)'},'fontname','helvetica');%,'fontweight','bold');
%     xlabel(hmei_list{j,1});
%     ylim([0 20])
%     
%     ix0 = (v2(:,2) <= prctile(v2(:,2), 25)); sum(ix0);
%     ix01 = (v2(:,2) > prctile(v2(:,2), 25)); sum(ix01);
%     [hh, pp_basal(j,1)] = ttest2(v1(ix0), v1(ix01));
%     
%     ix0 = (v2(:,2) <= prctile(v2(:,2), 50)); sum(ix0);
%     ix01 = (v2(:,2) > prctile(v2(:,2), 50)); sum(ix01);
%     [hh, pp_basal(j,2)] = ttest2(v1(ix0), v1(ix01));
end
r_basal= cell2table(r_basal); 
r_basal.Properties.RowNames= {hmei_list.cpd_name{1:8}};
% pp_basal= array2table(pp_basal);
% pp_basal.Properties.RowNames= {hmei_list.cpd_name{1:8}};
%disp(pp_basal) %t-test p-values for each drug

%% Last section uses grate_ccle_exp_soft_hdacsign, which was not calculated. Do not run this section.
%predicting variation in sensitivity between hme inhibitors based on basal 
%metabolic state and drug impact - comparison with drug sensitivity data from seashore-ludlow study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1 = grate_ccle_exp_soft_hdacsign(:, 5:8);
g1 = g1 - repmat(ignoreNaN(g1, @median,2),1,4); % 

[ix pos] = ismember(hcommon1, exptidcelllinemediamatch(:,1)); sum(ix) % match cell lines
v1 = hcommon_exptdat;
v1 = v1 - repmat(median(v1,2),1,4);
v2 = g1(pos,:);
v1 = v1(:);


ix0 = (v2(:) <= prctile(v2(:), 25)); sum(ix0)
ix01 = (v2(:) >= prctile(v2(:), 25)); sum(ix01)
[hh pp] = ttest2(v1(ix0), v1(ix01)) % 2e-52..

ix0 = (v2(:) <= prctile(v2(:), 50)); sum(ix0)
ix01 = (v2(:) >= prctile(v2(:), 50)); sum(ix01)
[hh pp] = ttest2(v1(ix0), v1(ix01)) % 3e-44

ix0 = (v2(:) <= 0.05); sum(ix0)
ix01 = (v2(:) >0.05); sum(ix01)
[hh pp] = ttest2(v1(ix0), v1(ix01)) %2e -62
mean(v1(ix0)) % -0.01
mean(v1(ix01)) % -2.61

    groups = NaN(size(ix0));
 groups(ix0) = 1; groups(ix01) = 2;
vv = NaN(length(v1), 2);

 vv(1:sum(groups == 1),1) = v1(groups == 1);
 vv(1:sum(groups == 2),2) = v1(groups == 2);
figure;
bh = boxplot(vv,'symbol','' ,'orientation','horizontal')
 set(gca,'yticklabel',{'No','Large'})
 xlabel('Observed differential sensitivity between drugs in a cell line')
title({'Predicted differential sensitivity between drugs' 'vs Observed differential sensitivity (AUC)'})


figure
h1 = histfit(vv(:,1),20,'kernel')%,
set(h1(1),'facecolor','g','facealpha',.15,'edgecolor','none')
set(h1(2),'color','g')
hold on
h2 = histfit(vv(:,2),20,'kernel')
set(h2(1),'facecolor','m','facealpha',.15,'edgecolor','none')
set(h2(2),'color','m')
hh = legend([h1(2),h2(2)],'No difference','Large difference')
title(hh,'Predicted differential acetylation')
