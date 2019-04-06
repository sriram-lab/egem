%% Methylation

% metabolic model with nuclear acetylation reaction
load supplementary_software_code acetylation_model 

% CCLE cell line names and gene expression data (z-transformed)
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz 

% cell line names, growth media, total bulk acetylation
load supplementary_software_code acetlevlistmedia acetlevellist acetlevellistval 

% methylation and acetylation levels of different sites from LeRoy et al.,
load methylation_proteomics_validation_data acet_meth_listval acet_meth_list_rowlab 

% list of nutrient conditions and uptake rates
load supplementary_software_code labels media_exchange1 mediareactions1

% 1 for maximizing SAM; 2 for maximizing methylated histones
meth_type = 1; 
MODE = 1;  % changed to rxn.
epsilon = 1E-2; rho = 1;
kappa = 1;
minfluxflag = 0; 

for i = 1:14
    iii = find(ismember(celllinenames_ccle1, acetlevellist(i)));
    if ~isempty(iii)
        iii  = iii(1);
        model2 = acetylation_model;
        
        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) > 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) < -2));
        
        % now set the media and glucose levels for different media conditions
        if ismember({'RPMI'} , acetlevlistmedia(i))
	        % no change rpmi
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5;
        elseif ismember({'DMEM'} , acetlevlistmedia(i))
	        % dmem
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;
        elseif ismember({'L15'} , acetlevlistmedia(i)) % NO GLUC AND LOW GAL
	        % L15   
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -0;
            % LOW GAL
            model2.lb(find(ismember(model2.rxns, {'EX_gal(e)'})))  = -0.9;
        elseif ismember({'McCoy 5A'} , acetlevlistmedia(i))
	        % mccoy
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*3/2;
        elseif ismember({'IMM'} , acetlevlistmedia(i))
	        % IMDM
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;
		end
        
		% Single gene deletion analysis
		[~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
		[~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
		disp(i)

		% Get the flux redistribution values associated with different media component addition and deletion
		[fluxstate_gurobi,grate_ccle_exp_acetdat(i,1), solverobj_ccle(i,1)] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[], minfluxflag);

		% Now let's add the methylation reaction we want; either maximize SAM production or histone methylation.
        if meth_type == 1
	        model2 = addReaction(model2,'amet_out',    'amet[n]  -> ');  % maximize nuclear SAM
	        rxnpos1  = [find(ismember(model2.rxns,'amet_out'));];
        elseif meth_type == 2
	        model2 = addReaction(model2,'EX_HistMET',    'Nmelys[n] -> '); % max nuclear methylation
            rxnpos1  = [find(ismember(model2.rxns,'EX_HistMET'));];
        end
		
		% limit methionine levels for all reactions in the model; it has to be non limiting
		[ix, pos]  = ismember({'EX_met_L(e)'},model2.rxns);
        model2.lb(pos) = -0.5; 
		model2.c(rxnpos1) = 1; % we're interested in this reaction
		
		% get the flux values
        [fluxstate_gurobi] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[],minfluxflag);
        grate_ccle_exp_acetdat(i,2) = fluxstate_gurobi(rxnpos1);
    end
end

% Calculate the correlation
[acetlevelcorr_amet, acetlevelcorrpv_amet ] = corr(grate_ccle_exp_acetdat(:,2), acet_meth_listval');
acetlevelcorr_amet = acetlevelcorr_amet';

% this has correlation with various methylation acetylation marks. 
figure;
barh([acetlevelcorr_amet], 1, 'edgecolor', 'w');
set(...
    gca, 'ytick', [1:length(acet_meth_list_rowlab)], ...
    'yticklabel', acet_meth_list_rowlab,...
    'fontsize', 8, ...
    'fontweight','bold');
set(gca,'TickDir', 'out');
set(gca,'box','off');
set(gca,'linewidth',2);
set(gcf,'color','white');
set(gca,'fontsize',12);
xlabel('Pearson Correlation');
ylabel('H3 methylation and acetylation positions');
xlim([-1,1]);
title('Correlation between histone mark expression and metabolic flux', 'fontweight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impact of growth media components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posgluc = 1385;  % glucose uptake reaction in recon1.
objpos = find(acetylation_model.c) %biomass objective
minfluxflag = 0; % no PFBA
epsilon_acetylation = 1; % higher weights for methylation compared to acetylation

for kappatype = 1:2
    if kappatype == 1, kappa  = 10; else kappa = 0.01;end
    
    for i = 1:50
        kappa1 = kappa;
        if (kappatype == 2) & (ismember(i,[2,3,5:19])) % trace elements
            kappa1 = kappa/100;
        elseif (kappatype == 1) & (ismember(i,[1;4])) % glucose or glutamine
            kappa1 = 3;
        end
        model2 = acetylation_model;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % change media..
        [ix pos]  = ismember({'EX_met_L(e)'},model2.rxns);
        model2.lb(pos) = -0.5; % it has to be non limiting...................
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if meth_type == 1
       % SAM flux maximize
                model2 = addReaction(model2,'amet_out',    'amet[n]  -> ');
                rxnpos  = find(ismember(model2.rxns,'amet_out'));
              elseif meth_type == 2
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model2 = addReaction(model2,'EX_HistMET',    'Nmelys[n] -> ');
        rxnpos  = find(ismember(model2.rxns,'EX_HistMET'));
              end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ix pos]  = ismember(mediareactions1(i), model2.rxns);
        model2.lb(pos) = -media_exchange1(i,1)*kappa1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [solf.x,sol11] =  constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);
        
        str = ['media_change_growth_',num2str(kappatype),'(i,1) = solf.x(objpos);'];
        if ~isempty(solf.x) &  ~isnan(solf.x)
            eval(str)
        end
        
        j = 1;
        model3 = model2;
        model3.c(rxnpos) = epsilon_acetylation;
        [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
        str = ['media_change_histone_acet_nuc_',num2str(kappatype),'(i,j) = solf.x(rxnpos);'];
        if ~isempty(solf.x) &  ~isnan(solf.x)
            eval(str)
        end
        disp(i)
    end
    
    
    disp(kappatype)
end

labels(2) = {'Glutathione'};
idx = [1:4,11:50];
figure;
bar([media_change_histone_acet_nuc_1(idx,1) media_change_histone_acet_nuc_2(idx,1) ],1,'edgecolor','w');
title('Methylation levels in different growth conditions','fontweight','bold')
set(gca,'xtick',[1:length(mediareactions1(idx))],'xticklabel',labels(idx),'fontsize',8,'fontweight','bold','XTickLabelRotation',45)
set(gca,'TickDir', 'out')
set(gca,'box','off')
set(gca,'linewidth',2)
set(gcf,'color','white')
set(gca,'fontsize',12)
%ylabel('SAM- Flux') % if Maximizing SAM
ylabel('Methyl- Flux')
h = legend({'Excess','Depletion'});



%% Methylation and methylation levels 





