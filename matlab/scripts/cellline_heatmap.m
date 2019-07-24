for i = 1:14
    iii = find(ismember(ccle_names, leroy_names));
    if ~isempty(iii)
        iii = iii(1);
        model2 = model;
        
        %find up and down-regulated genes in each cell line
        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) >= 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) <= -2));
        
        % Set the glucose uptake based on media; Default glucose is -5 for rpmi
        if ismember({'RPMI'} , cellmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5;% no change rpmi
        elseif ismember({'DMEM'} , cellmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;% dmem.. 
        elseif ismember({'L15'} , cellmedia(i))   % NO glucose.. LOW Galactose
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -0;% L15
            model2.lb(find(ismember(model2.rxns, {'EX_gal(e)'})))  = -0.9;%
        elseif ismember({'McCoy 5A'} , cellmedia(i)) 
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*3/2;% mccoy
        elseif ismember({'IMM'} , cellmedia(i))
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'})))  = -5*4.5/2;% IMDM
        end
        
        % Find reactions from differentially expressed genes
        [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
        [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
        
        disp([i,' media conditions done'])
        
        [fluxstate_gurobi, ccle_grate(i,1), solverobj_ccle(i,1)] =...
            constrain_flux_regulation(model2, onreactions, offreactions,...
            kappa, rho, epsilon, MODE, [], minfluxflag);
        model2.c(rxnpos) = epsilon2;
        
        % Recalc the flux value again when we set the objective coefficient
        % for the rxn of interest to epsilon2. Get the flux value. 
        [fluxstate_gurobi] =  constrain_flux_regulation(model2, onreactions,...
            offreactions, kappa, rho, epsilon, MODE,[], minfluxflag);
        ccle_grate(i,2) = fluxstate_gurobi(rxnpos);
    end
end

for x = 1:length(labels)
    [R_up, pval_up] = corr(ccle_grate(x,2), media_xchange_1(:, 1));
    [R_down, pval_down] = corr(ccle_grate(x,2), media_xchange_2(:, 1));