      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % this code reproduces all the main text figures from the manuscript (Shen et al)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% impact of nutrient sources on acetylation - figure 2A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load recon1
%load supplementary_software_code acetylation_model %contains metabolic model with nuclear acetylation reaction
load supplementary_software_code labels media_exchange1 mediareactions1 %list of nutrient conditions and uptake rates
posgluc = 1385;  % glucose uptake reaction in recon1. 
% LF changed reaction to find
 rxnpos  = find(ismember(metabolicmodel.rxns,'LYSMTF1n')); % rxnpos = 2451
objpos = find(metabolicmodel.c); %biomass objective
minfluxflag = 0; % no PFBA
%epsilon_acetylation = 1E-3;
epsilon_methylation = 1E-2; % or 1E-1

    for kappatype = 1:2
        if kappatype == 1, kappa  = 10; else kappa = 0.01; end 
        %kappatype=1 means high  [nutrient]. kappytype=2 means low
        for i = 1:50
            kappa1 = kappa;
            if (kappatype == 2) && (ismember(i,[2,3,5:19])) % trace elements
                kappa1 = kappa/100;
            elseif (kappatype == 1) && (ismember(i,[1;4])) % glucose or glutamine
                kappa1 = 3;
            end
            model2 = metabolicmodel;
            % change media..
            [ix, pos] = ismember(mediareactions1(i), model2.rxns);
            model2.lb(pos) = -media_exchange1(i,1)*kappa1;
            
            [solf.x,sol11] = constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);
            
            str = ['media_change_growth_',num2str(kappatype),'(i,1) = solf.x(objpos);'];
            if ~isempty(solf.x) && all(~isnan(solf.x))
                eval(str)
            end
            
            j = 1;
            model3 = model2;
            model3.c(rxnpos) = epsilon_methylation;
            [solf.x,sol11] = constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
            str = ['media_change_histone_acet_nuc_',num2str(kappatype),'(i,j) = solf.x(rxnpos);'];
            if ~isempty(solf.x) && all(~isnan(solf.x))
                eval(str)
            end
            disp(i)
        end
        
        
        disp(kappatype)
    end
    
    labels(2) = {'Glutathione'};
    idx = [1:4,20:50];
    figure;
    bar([media_change_histone_acet_nuc_1(idx,1) media_change_histone_acet_nuc_2(idx,1) ],1,'edgecolor','w');
    title('Acetylation levels in different growth conditions','fontweight','bold')
    set(gca,'xtick',[1:length(mediareactions1(idx))],'xticklabel',labels(idx),'fontsize',8,'fontweight','bold','XTickLabelRotation',45)
    set(gca,'TickDir', 'out')
    set(gca,'box','off')
    set(gca,'linewidth',2)
    set(gcf,'color','white')
    set(gca,'fontsize',12)
    ylabel('Acetyl- Flux')
h = legend({'Excess','Depletion'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% impact of gene deletion on acetylation - figure 2B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            unqgenes = unique(acetylation_model.genes); %all genes in the model

    for i = 2:length(unqgenes)
        
        modeltemp = acetylation_model ;
         modeltemp = deleteModelGenes(modeltemp,unqgenes(i));
        model3 = modeltemp;
      
        
        [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
        if ~isempty(solf.x) && ~isnan(solf.x)
            acet_screen_rpmi(i,2) = solf.x(objpos); % impact on growth
        end
        
        model3.c(rxnpos) = epsilon_acetylation; % epsilon is a weight
        [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
        if ~isempty(solf.x) && ~isnan(solf.x)
            acet_screen_rpmi(i,1) = solf.x(rxnpos); % impact on acetylation flux
        end
         
        disp(i)
    end
    
    acet_screen_rpmi(1,:) = NaN;

     [sx spos] = sort(acet_screen_rpmi(:,1));
 acetgenessorted = unqgenes(spos);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% impact of culture media components on acetylation - figure 2C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glucflag = logical([1 0 1 1 1 1 0 0 0 0]);
aaflag =  logical([1 0 1 0 1 0 0 1 0 1]);
pyrflag = logical([1 0 1 1 0 0 0 0 1 1]);
glnflag = logical([1 0 0 1 1 0 1 0 0 1]);
expval = [1 0 1 1 1 1 1 0 1 1];

% part 2
glucflag1 = [ ones(1,6),  0 ,0, 0,  0];
expval1 = [ ones(1,7),  0 , 1,  0];


[ix, aapos] = ismember( mediareactions1(20:37),acetylation_model.rxns);
posgluc = 1385;  % glucose uptake reaction in recon1.
glnpos = 1386; % glutamine
pyrpos = 1509; % pyruvate
acetatepos = 1238; % acetate
fattyacidpos = 1445; % linoleic acid

model3 = acetylation_model;
model3.c(rxnpos) = epsilon_acetylation ;
[solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
wild_type_acet = solf.x(rxnpos); % default acetylation flux

for i = 1:10
    model3 = acetylation_model;
    if ~glucflag(i)
        model3.lb(posgluc) = 0;
    end
    
    if ~aaflag(i)
        model3.lb(aapos) = 0;
    end
    
    if ~glnflag(i)
        model3.lb(glnpos) = 0;
    end
    
    if pyrflag(i)
        model3.lb(pyrpos) = -10;
    end
    
    
    model3.c(rxnpos) = epsilon_acetylation;
    [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    acet_media_screen_dmem(i,1) = solf.x(rxnpos); % impact on acetylation flux
    

end
disp('consistency with 10 different experimental conditions - Part I')
sum((acet_media_screen_dmem > 0.05*wild_type_acet) == expval') % 9



for i = 1:10
    model3 = acetylation_model;
    if ~glucflag1(i)
        model3.lb(posgluc) = 0;
    end
   
    
    switch i 
        case 2 % no ca+
                    model3.lb(ismember( {'EX_ca2(e)'},acetylation_model.rxns)) = 0;
        case 4 % no Phosphate
                    model3.lb(ismember( {'EX_pi(e)'},acetylation_model.rxns)) = 0;
        case 6 % no vitamins in dmem.. compmosition from sigma 
                    model3.lb(ismember( {'EX_thm(e)';'EX_ribflv(e)';'EX_pydx(e)';'EX_ncam(e)';'EX_inost(e)';'EX_5mthf(e)';'EX_pnto_R(e)';'EX_chol(e)'},acetylation_model.rxns)) = 0;
        case 7 %  add acetate
                    model3.lb(acetatepos) = -20;
        case 9 % add fatty acid
                model3.lb(fattyacidpos) = -5;
    end

        
    model3.c(rxnpos) = epsilon_acetylation;
    [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    acet_media_screen_dmem1(i,1) = solf.x(rxnpos); % impact on acetylation flux
end

disp('consistency with 10 different experimental conditions - Part II')
sum((acet_media_screen_dmem1 > 0.05*wild_type_acet) == expval1') % 10

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% impact of basal metabolic state of cell lines on acetylation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load supplementary_software_code celllinenames_ccle1   ccleids_met ccle_expression_metz %contains CCLE cell line names, gene expression data (z-transformed)
load supplementary_software_code acetlevlistmedia acetlevellist acetlevellistval %contains cell line names, growth media , total bulk acetylation
MODE = 1;  % reaction (1) or gene list (0)
epsilon = 1E-2; rho = 1;
kappa = 1; 
minfluxflag = 0; % no PFBA

 for i = 1:14
    % match cell line in CCLE data
    iii = find(ismember(celllinenames_ccle1, acetlevellist(i)));
    if ~isempty(iii)
        iii  = iii(1);
        model2 = acetylation_model;
                 %find up and down-regulated genes in each cell line
        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) > 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) < -2));
        
        % set the glucose uptake based on media
        % default glucose is -5 for rpmi
        if ismember({'RPMI'} , acetlevlistmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5;% no change rpmi
        elseif ismember({'DMEM'} , acetlevlistmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;% dmem.. 
        elseif ismember({'L15'} , acetlevlistmedia(i))   % NO glucose.. LOW Galactose
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -0;% L15
            model2.lb(find(ismember(model2.rxns, {'EX_gal(e)'})))  = -0.9;%
        elseif ismember({'McCoy 5A'} , acetlevlistmedia(i)) 
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*3/2;% mccoy
        elseif ismember({'IMM'} , acetlevlistmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;% IMDM
       
        end
        
              %find reactions from differentially expressed genes
        [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
        [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
        
        disp(i)
       
            [fluxstate_gurobi,grate_ccle_exp_acetdat(i,1), solverobj_ccle(i,1)] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[], minfluxflag);  % impact on growth
            
            model2.c(rxnpos) = epsilon_acetylation;
            [fluxstate_gurobi] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[],minfluxflag);
            grate_ccle_exp_acetdat(i,2) = fluxstate_gurobi(rxnpos); %acetylation flux
            
    end
 end
 
     figure;
    plot(grate_ccle_exp_acetdat(:,2), acetlevellistval(1,:)','o','markerfacecolor',[ 0.9020    0.3804    0.0039],'markeredgecolor','k') 
fb = polyfit(grate_ccle_exp_acetdat(:,2), acetlevellistval(1,:)',1);
fb1 = polyval(fb,grate_ccle_exp_acetdat(:,2));
hold on; plot(grate_ccle_exp_acetdat(:,2),fb1,'r-','linewidth',2);
text(grate_ccle_exp_acetdat(:,2) + 0.2,acetlevellistval(1,:)',acetlevellist)%,'horizontalalignment','right','fontsize',10,'fontname','helvetica','fontweight','bold')
xlabel('Acetylation flux');ylabel('Bulk H3K9 Acetylation')

  [acetlevelcorr acetlevelcorrpv ] = corr(grate_ccle_exp_acetdat(:,2), acetlevellistval(1,:)') %correlation with h3k9 acetylation
 % [acetlevelcorr acetlevelcorrpv ] = corr(grate_ccle_exp_acetdat(:,2), acetlevellistval(1,:)','type','spearman') %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% impact of nutrient environment of cell lines on sensitivity to deacetylase inhibitor - vorinostat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 load supplementary_software_code recon1biologpm1match biologratio biolognewpm196 %contains BIOLOG phenotype array data

for i = 1:96
      model2 = acetylation_model;
    if ~isnan(recon1biologpm1match(i))
        model2.lb(posgluc) = -0.1;
        model2.lb(recon1biologpm1match(i)) = -10;
        
        [solf.x,sol11] =  constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);
        disp(i)
        
        growth_basal2(i,1) = solf.x(objpos);  % impact on growth
        model2.c(rxnpos) = epsilon_acetylation;
        [solf.x,sol11] =  constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);
        growth_basal2(i,2) = solf.x(rxnpos); %acetylation flux
        
    else
        growth_basal2(i,:) = NaN;
    end
end

ixxn2 = ~isnan( growth_basal2(:,1)) & biolognewpm196(:,1) > 2; sum(ixxn2) % 
xx = biolognewpm196(ixxn2,2)./biolognewpm196(ixxn2,1);
tel = [2:3,5,7:19];
yy = [growth_basal2(ixxn2,:)];
[hh pp] = corr([yy(tel,1) xx(tel) ],'type','spearman') % . -0.36 %correlation between predicted growth and vorinostat sensitivity
[hh pp] = corr([yy(tel,2) xx(tel) ],'type','spearman') % -0.67 %correlation between predicted acetylation flux and vorinostat sensitivity


figure; plot(xx(tel), yy(tel,2) ,'o','markerfacecolor',[ 0.9020    0.3804    0.0039],'markeredgecolor','k') % 0.5
xlim([0.5 2])
fb = polyfit(xx(tel), yy(tel,2),1);
fb1 = polyval(fb,xx(tel));
hold on; plot(xx(tel),fb1,'r-','linewidth',2);
ylabel('Acetylation flux');xlabel('Vorinostat treatment vs control ratio from Biolog arrays')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% impact of basal metabolic state of CCLE cell lines on sensitivity to demethylase inhibitors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load supplementary_software_code  ctd2celllineidname exptidcelllinemediamatch ctd2celllineidname_id r* 
%data from seashore-ludlow study, contains cell line names, growth media
load supplementary_software_code  ctd2compoundidname_id drug_auc_expt ctd2compoundidname_name 
%data from seashore-ludlow study, contains drug names, drug sensitivity data
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz  
%contains CCLE cell line names, gene expression data (z-transformed)
load supplementary_software_code hcommon_exptdat hcommon1 
%data from seashore-ludlow study, contains drug names, drug sensitivity data for cell lines that were screened against all 4 hdac inhibitors
load supplementary_software_code hdacexpfcs hdacexpallgeneids 
%contains gene expression data after treatment with hdac inhibitors
rxnpos = 2451; % loading the above variables changes rxnpos to 3754, which is out of bounds for fluxstate_gurobi

MODE = 1;  % reaction (1) or gene list (0)
epsilon = 1E-2; rho = 1;
kappa = 1; % parameters for integrating transcriptomics data. kappa is the strength of down regulation of genes (Chandrasekaran & Price, PNAS, 2010)
minfluxflag = 0; % no Pfba 
hdactransint = 1;
basalflag = 1;

grate_ccle_exp_soft = NaN(length(exptidcelllinemediamatch),2); %contains acetylation flux based on basal metabolic state
grate_ccle_exp_soft_hdacsign = NaN(length(exptidcelllinemediamatch),8); %contains acetylation flux based on basal metabolic state and impact of each hdac inhibitor teatment

for i = 1:length(exptidcelllinemediamatch)
    % match cell line data in CCLE with CTD2
    ii = find(ismember(ctd2celllineidname_id,  exptidcelllinemediamatch(i,2)));
    iii = find(ismember(celllinenames_ccle1, ctd2celllineidname(ii,1)));
    if ~isempty(iii)
        iii  = iii(1);
         model2 = metabolicmodel;
      
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
        [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
        [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
%        
        
        disp(i)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get basal metabolic state based on transcriptome
            if basalflag
        [fluxstate_gurobi,grate_ccle_exp_soft(i,1), solverobj_ccle(i,1)] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[], minfluxflag); % impact on growth
        model2.c(rxnpos) = epsilon_methylation ;
        [fluxstate_gurobi] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE,[],minfluxflag);
        grate_ccle_exp_soft(i,2) = fluxstate_gurobi(rxnpos); %acetylation flux
           end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % hdac inhibitor impact on transcriptome
        if hdactransint
            for kk = 4:-1:1   % match drug order
                ongenesh = hdacexpallgeneids((hdacexpfcs(:,kk) > 1.1) & (hdacexpfcs(:,kk + 4) < 0.01));
                ongenesh = intersect(ongenesh, model2.genes);
                offgenesh = hdacexpallgeneids((hdacexpfcs(:,kk) < 0.9) & (hdacexpfcs(:,kk + 4) < 0.01));
                offgenesh = intersect(offgenesh, model2.genes);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,~,onreactionsh,~] =  deleteModelGenes(model2, ongenesh);
                [~,~,offreactionsh,~] =  deleteModelGenes(model2, offgenesh);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % remove conflicting transcripts
                onreactionsh0 = setdiff(onreactionsh, offreactions);
                offreactionsh0 = setdiff(offreactionsh, onreactions);
                onreactions0 = setdiff(onreactions, offreactionsh);
                offreactions0 = setdiff(offreactions, onreactionsh);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                onreactions1h = [onreactions0;  onreactionsh0];
                offreactions1h = [offreactions0;  offreactionsh0];
                epsilon2 = [repmat(0, size(offreactions1h))];
                model2.c(rxnpos) = 0;
         %        [fluxstate_gurobi,grate_ccle_exp_soft_hdacsign(i,kk)] =  constrain_flux_regulation(model2,onreactions1h,offreactions1h,kappa,rho,epsilon,MODE, epsilon2,minfluxflag); % impact on growth
                model2.c(rxnpos) = epsilon_methylation ;
                  [fluxstate_gurobi] =  constrain_flux_regulation(model2,onreactions1h,offreactions1h,kappa,rho,epsilon,MODE,epsilon2,minfluxflag);
                 grate_ccle_exp_soft_hdacsign(i,kk + 4) = fluxstate_gurobi(rxnpos); %acetylation flux
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end

    else
        grate_ccle_exp_soft(i,:) = NaN;
        grate_ccle_exp_soft_hdacsign(i,:) = NaN;
    end
end


figure; h = histogram(grate_ccle_exp_soft(:,2),70);
 xlabel('Predicted acetylation flux')
 ylabel('Total cell lines')
title('Distribution of acetylation flux among CCLE cell lines','fontweight','normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%predicting sensitivity to hdac inhibitors based on basal metabolic state - comparison with drug sensitivity data from seashore-ludlow study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdaclist = {'LBH-589','vorinostat','entinostat','belinostat'}
%hdmelist =

for j = 1:4
fx = find(ismember(ctd2compoundidname_name, hdaclist(j))) 
ix = ismember(drug_auc_expt(:,3), ctd2compoundidname_id(fx)); 
sum(ix) % 847. ix is a 1D array. Each drug appears multiple times (appears for each medium in drug_auc_expt)
hdac_auc_dat = drug_auc_expt(ix,:);

[ix pos] = ismember(hdac_auc_dat(:,1), exptidcelllinemediamatch(:,1)); sum(ix) %
v1 = hdac_auc_dat(:,2);
v2 = grate_ccle_exp_soft(pos,:);

ix0 = (v2(:,2) < 0.05);sum(ix0) % 0.05 is  the default threshold value 
ix01 = (v2(:,2) > 0.05);sum(ix01)
[hh pp_basal(j,3)] = ttest2(v1(ix0), v1(ix01)); %

  groups = NaN(size(ix0));
 groups(ix0) = 1; groups(ix01) = 2;
 vv = NaN(length(v1), 2);
 vv(1:sum(groups == 1),1) = v1(groups == 1);
 vv(1:sum(groups == 2),2) = v1(groups == 2);

 figure;
 %clf; UnivarScatter(vv,'Width', 0.3, 'PointSize', 11,'MarkerEdgeColor','w','LineWidth',0.1);%, 'markerfacealpha',0.5);
 %hold on; 
 bh = boxplot(v1, groups,'symbol','')
 set(gca,'xticklabel',{'Low flux','High flux'})
ylabel({'Sensitivity (AUC)'},'fontname','helvetica');%,'fontweight','bold');
xlabel(hdaclist(j));
ylim([0 20])

ix0 = (v2(:,2) <= prctile(v2(:,2), 25)); sum(ix0)
ix01 = (v2(:,2) > prctile(v2(:,2), 25)); sum(ix01)
[hh pp_basal(j,1)] = ttest2(v1(ix0), v1(ix01)) ;% 

ix0 = (v2(:,2) <= prctile(v2(:,2), 50)); sum(ix0)
ix01 = (v2(:,2) > prctile(v2(:,2), 50)); sum(ix01)
[hh pp_basal(j,2)] = ttest2(v1(ix0), v1(ix01)) ;% 

end
disp(pp_basal) %t-test p-values for each drug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%predicting variation in sensitivity between hdac inhibitors based on basal 
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
