%% Methylation
% Arguments for maximizing flux for different metabolites (demand reactions):
    % 1: SAM
    % 2: Histone Methylation 1
    % 3: Histone Methylation 2
    % 4: Histone Methylation 3
    % 5: Folate
    % 6: Choline
    % 7: DNA Methylation
    % 8: Serine
    % 9: Glycine
    % 10: THF
    % 11: AKG
    % 12: Succ
    % 13: Fum
    % 14: ATP
    % 15: NAD
    % 16: NADH
    % 17: ADP

% Arguments for compartments
    % c: cyto
    % m: mito
    % n: nuc

%model = model;
%model = acetylation_model;
%meth_type = 16;

% Finding specific reactions in the model. Otherwise, comment it out.
rxnpos1  = [find(ismember(model.rxns, 'MAT2n'));];
nam = 'MAT2n';

%compartment = 'n';
MODE = 1;  % changed to rxn.
epsilon = 1; 
rho = 1;
kappa = 1;
minfluxflag = 0; 

%% See if there is a correlation between gene expression and flux
for i = 1:14
    iii = find(ismember(celllinenames_ccle1, acetlevellist(i)));
    if ~isempty(iii)
        iii  = iii(1);
        model2 = model;

        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) >= 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) <= -2));

        % now set the media and glucose levels for different media conditions
        if ismember({'RPMI'} , acetlevlistmedia(i))
            % no change rpmi
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5;
        elseif ismember({'DMEM'} , acetlevlistmedia(i))
            % dmem
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*4.5/2;
        elseif ismember({'L15'} , acetlevlistmedia(i)) % NO GLUC AND LOW GAL
            % L15   
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -0;
            % LOW GAL
            model2.lb(find(ismember(model2.rxns, {'EX_gal(e)'}))) = -0.9;
        elseif ismember({'McCoy 5A'} , acetlevlistmedia(i))
            % mccoy
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*3/2;
        elseif ismember({'IMM'} , acetlevlistmedia(i))
            % IMDM
            model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*4.5/2;
        end

        % Single gene deletion analysis
        [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
        [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
        disp(i)

        % Get the flux redistribution values associated with different media component addition and deletion
        [fluxstate_gurobi, grate_ccle_exp_dat(i,1),  solverobj_ccle(i,1)] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[], minfluxflag);

        % Now let's add the demand reaction we want
        if (~exist('meth_type','var')) || (isempty(meth_type))
            rxnpos1 = rxnpos1;
        elseif meth_type == 1
            model2 = addReaction(model2, 'DM_amet', 'reactionFormula', ['amet[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_amet'));];
            nam = 'DM_amet';
        elseif meth_type == 2
            model2 = addReaction(model2, 'EX_HistMET1', 'reactionFormula', ['Nmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'EX_HistMET1'));];
            nam = 'EX_HistMET1';
        elseif meth_type == 3
            model2 = addReaction(model2, 'DM_HistMET2', 'reactionFormula', ['Ndmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_HistMET2'));];
            nam = 'DM_HistMET2';
        elseif meth_type == 4
            model2 = addReaction(model2, 'DM_HistMET3', 'reactionFormula', ['Ntmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_HistMET3'));];
            nam = 'DM_HistMET3';
        elseif meth_type == 5
            model2 = addReaction(model2, 'DM_fol', 'reactionFormula', ['fol[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_fol'));];
            nam = 'DM_fol';
        elseif meth_type == 6
            model2 = addReaction(model2, 'DM_chol', 'reactionFormula', ['chol[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_chol'));];
            nam = 'DM_chol';
        elseif meth_type == 7
            model2 = addReaction(model2, 'DM_DNAMe', 'reactionFormula', ['dna5mtc[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_DNAMe'));];
            nam = 'DM_DNAMe';
        elseif meth_type == 8
            model2 = addReaction(model2, 'DM_ser', 'reactionFormula', ['ser[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_ser'));];
            nam = 'DM_ser';
        elseif meth_type == 9
            model2 = addReaction(model2, 'DM_gly', 'reactionFormula', ['gly[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_gly'));];
            nam = 'DM_gly';
        elseif meth_type == 10
            model2 = addReaction(model2, 'DM_thf', 'reactionFormula', ['thf[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_thf'));];
            nam = 'DM_thf';
        elseif meth_type == 11
            model2 = addReaction(model2, 'DM_akg', 'reactionFormula', ['akg[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_akg'));];
            nam = 'DM_akg';
        elseif meth_type == 12
            model2 = addReaction(model2, 'DM_succ', 'reactionFormula', ['succ[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_succ'));];
            nam = 'DM_succ';
        elseif meth_type == 13
            model2 = addReaction(model2, 'DM_fum', 'reactionFormula', ['fum[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_fum'));];
            nam = 'DM_fum';
        elseif meth_type == 14
            model2 = addReaction(model2, 'DM_atp', 'reactionFormula', ['atp[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_atp'));];
            nam = 'DM_atp';
        elseif meth_type == 15
            model2 = addReaction(model2, 'DM_nad', 'reactionFormula', ['nad[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_nad'));];
            nam = 'DM_nad';
        elseif meth_type == 16
            model2 = addReaction(model2, 'DM_nadh', 'reactionFormula', ['nadh[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_nadh'));];
            nam = 'DM_nadh';
        elseif meth_type == 17
            model2 = addReaction(model2, 'DM_adp', 'reactionFormula', ['adp[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_adp'));];
            nam = 'DM_adp';
        end

        % limit methionine levels for all reactions in the model; it has to be non limiting
        [ix, pos]  = ismember({'EX_met_L(e)'}, model2.rxns);
        model2.lb(pos) = -0.5;
        model2.c(3743) = 0;
        model2.c(rxnpos1) = 0.01; % we're interested in this reaction

        % get the flux values from iMAT
        [fluxstate_gurobi] =  constrain_flux_regulation(model2,...
            onreactions, offreactions,...
            kappa, rho, epsilon,...
            MODE ,[], minfluxflag);
        grate_ccle_exp_dat(i,2) = fluxstate_gurobi(rxnpos1);
    end
end

% Calculate the correlation
[acetlevelcorr_amet, acetlevelcorrpv_amet ] = corr(grate_ccle_exp_dat(:,2), acet_meth_listval');
acetlevelcorr_amet = acetlevelcorr_amet';

% Make plot of correlation
fig = figure;
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
set(gcf, 'Position', [100, 100, 700, 800])
xlabel('Pearson Correlation');
ylabel('H3 methylation and acetylation positions');
xlim([-1,1]);
title(['Correlation between histone mark expression and ' nam 'metabolic flux'], 'fontweight', 'bold');
saveas(fig(1), ['./../figures/fig/ac-model' nam '-corr-memodel.fig']);
saveas(fig(1), ['./../figures/tiff/ac-model' nam '-corr-memodel.tif']);

%% Determine the impact of excess or depleting growth media components on flux

posgluc = 1385;  % glucose uptake reaction in RECON1
objpos = find(model.c) %biomass objective
minfluxflag = 0; % no PFBA
new_epsilon = 1; % higher weights for methylation compared to acetylation

for kappatype = 1:2
    if kappatype == 1, kappa  = 10; else kappa = 0.01;end

    for i = 1:50
        kappa1 = kappa;
        if (kappatype == 2) & (ismember(i,[2,3,5:19])) % trace elements
            kappa1 = kappa/100;
        elseif (kappatype == 1) & (ismember(i,[1;4])) % glucose or glutamine
            kappa1 = 3;
        end
        model2 = model;

        % change media
        [ix, pos]  = ismember({'EX_met_L(e)'},model2.rxns);
        model2.lb(pos) = -0.5; % it has to be non limiting
        
        if (~exist('meth_type','var')) || (isempty(meth_type))
            rxnpos1 = rxnpos1;
        elseif meth_type == 1
            model2 = addReaction(model2, 'DM_amet', 'reactionFormula', ['amet[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_amet'));];
            nam = 'DM_amet';
        elseif meth_type == 2
            model2 = addReaction(model2, 'EX_HistMET1', 'reactionFormula', ['Nmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'EX_HistMET1'));];
            nam = 'EX_HistMET1';
        elseif meth_type == 3
            model2 = addReaction(model2, 'DM_HistMET2', 'reactionFormula', ['Ndmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_HistMET2'));];
            nam = 'DM_HistMET2';
        elseif meth_type == 4
            model2 = addReaction(model2, 'DM_HistMET3', 'reactionFormula', ['Ntmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_HistMET3'));];
            nam = 'DM_HistMET3';
        elseif meth_type == 5
            model2 = addReaction(model2, 'DM_fol', 'reactionFormula', ['fol[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_fol'));];
            nam = 'DM_fol';
        elseif meth_type == 6
            model2 = addReaction(model2, 'DM_chol', 'reactionFormula', ['chol[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_chol'));];
            nam = 'DM_chol';
        elseif meth_type == 7
            model2 = addReaction(model2, 'DM_DNAMe', 'reactionFormula', ['dna5mtc[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_DNAMe'));];
            nam = 'DM_DNAMe';
        elseif meth_type == 8
            model2 = addReaction(model2, 'DM_ser', 'reactionFormula', ['ser[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_ser'));];
            nam = 'DM_ser';
        elseif meth_type == 9
            model2 = addReaction(model2, 'DM_gly', 'reactionFormula', ['gly[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_gly'));];
            nam = 'DM_gly';
        elseif meth_type == 10
            model2 = addReaction(model2, 'DM_thf', 'reactionFormula', ['thf[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_thf'));];
            nam = 'DM_thf';
        elseif meth_type == 11
            model2 = addReaction(model2, 'DM_akg', 'reactionFormula', ['akg[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_akg'));];
            nam = 'DM_akg';
        elseif meth_type == 12
            model2 = addReaction(model2, 'DM_succ', 'reactionFormula', ['succ[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_succ'));];
            nam = 'DM_succ';
        elseif meth_type == 13
            model2 = addReaction(model2, 'DM_fum', 'reactionFormula', ['fum[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_fum'));];
            nam = 'DM_fum';
        elseif meth_type == 14
            model2 = addReaction(model2, 'DM_atp', 'reactionFormula', ['atp[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_atp'));];
            nam = 'DM_atp';
        elseif meth_type == 15
            model2 = addReaction(model2, 'DM_nad', 'reactionFormula', ['nad[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_nad'));];
            nam = 'DM_nad';
        elseif meth_type == 16
            model2 = addReaction(model2, 'DM_nadh', 'reactionFormula', ['nadh[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_nadh'));];
            nam = 'DM_nadh';
        elseif meth_type == 17
            model2 = addReaction(model2, 'DM_adp', 'reactionFormula', ['adp[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_adp'));];
            nam = 'DM_adp';
        end

        [ix, pos]  = ismember(mediareactions1(i), model2.rxns);
        model2.lb(pos) = -media_exchange1(i,1)*kappa1;

        [solf.x, sol11] =  constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);

        str = ['media_change_growth_',num2str(kappatype),'(i,1) = solf.x(objpos);'];
        if ~isempty(solf.x) & ~isnan(solf.x)
            eval(str)
        end

        j = 1;
        model3 = model2;
        model3.c(rxnpos1) = new_epsilon;
        [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
        str = ['media_change_histone_acet_nuc_',num2str(kappatype),'(i,j) = solf.x(rxnpos1);'];
        if ~isempty(solf.x) &  ~isnan(solf.x)
            eval(str)
        end
        disp(i)
    end
    disp(kappatype)
end

labels(2) = {'Glutathione'};
idx = [1:4,20:50];

fig = figure;
data = [media_change_histone_acet_nuc_1(idx, 1), media_change_histone_acet_nuc_2(idx, 1)];
plt = barh(data, 'edgecolor', 'w');
set(plt(2), 'FaceColor', hex2rgb('#C6393D'));
set(plt(1), 'FaceColor', hex2rgb('#BDCD6C'));
title('Varying media components', 'fontweight', 'bold');
set(gca,'ytick', [1:length(mediareactions1(idx))], 'yticklabel',...
    labels(idx), 'fontsize', 8, 'fontweight', 'bold');
set(gca,'TickDir', 'out');
set(gca,'box', 'off');
set(gca,'linewidth', 2);
set(gcf,'color', 'white');
set(gca,'fontsize', 12);
set(gcf, 'Position', [100, 100, 700, 800])
xlabel([nam ' flux (mmol/gDW*hr)']);
h = legend({'Excess', 'Depletion'});
legend boxoff;
saveas(fig(1), ['./../figures/fig/ac-model' nam '-media-memodel.fig']);
saveas(fig(1), ['./../figures/tiff/ac-model' nam '-media-memodel.tif']);






