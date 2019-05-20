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
%model = metabolicmodel;
%model = acetylation_model;

model_nam = 'methyl';

% Finding specific reactions in the model. Otherwise, comment it out.
%rxnpos1  = [find(ismember(model.rxns, 'DM_HistMe'));];
nam = 'SAM';
meth_type = 1;
compartment = 'n';
MODE = 1;  % changed to rxn.
epsilon = 1E-2; 
rho = 1;
kappa = 1E-3;
minfluxflag = 0; 

%% See if there is a correlation between gene expression and flux
for i = 1:14
    iii = find(ismember(celllinenames_ccle1, acetlevellist(i)));
    if ~isempty(iii)
        iii  = iii(1);
        model2 = model;

        ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) >= 2));
        offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) <= -2));

        % Set the media and glucose levels for different media conditions
        % corresponding to each cancer cell line.
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
        [fluxstate_gurobi, grate_ccle_exp_dat(i,1), solverobj_ccle(i,1)] =...
            constrain_flux_regulation(model2, onreactions, offreactions,...
            kappa, rho, epsilon, MODE, [], minfluxflag);

        % Now let's add the demand reaction we want
        if (~exist('meth_type','var')) || (isempty(meth_type))
            rxnpos1 = rxnpos1;
        elseif meth_type == 1
            model2 = addReaction(model2, 'DM_amet', 'reactionFormula', ['amet[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_amet'));];
            nam = 'DM_amet';
        elseif meth_type == 2
            model2 = addReaction(model2, 'DM_HistMET1', 'reactionFormula', ['Nmelys[' compartment '] -> ']);
            rxnpos1  = [find(ismember(model2.rxns, 'DM_HistMET1'));];
            nam = 'DM_HistMET1';
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
            model2 = addReaction(model2, 'DM_ser', 'reactionFormula', ['ser-L[' compartment '] -> ']);
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
        model2.c(3743) = 0; % Force the biomass function to be 0
        model2.c(rxnpos1) = 0.01; % Set the objective coefficient to be low.

        % Get the constrained flux values from iMAT
        [fluxstate_gurobi] =  constrain_flux_regulation(model2,...
            onreactions, offreactions, kappa, rho, epsilon,...
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
%title(['Correlation between histone mark expression and ' nam ' metabolic flux'], 'fontweight', 'bold');
saveas(fig(1), ['./../figures/fig/' model_nam '-' nam '-corr.fig']);
saveas(fig(1), ['./../figures/tiff/' model_nam '-' nam '-corr.tif']);

%% Determine the impact of excess or depleting growth media components on flux

posgluc = 1385;  % glucose uptake reaction in RECON1
objpos = find(model.c); % biomass objective
minfluxflag = 0; 
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
            nam = 'SAM consumption';
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
            model2 = addReaction(model2, 'DM_ser', 'reactionFormula', ['ser-L[' compartment '] -> ']);
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
%title('Varying media components', 'fontweight', 'bold');
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
saveas(fig(1), ['./../figures/fig/' model_nam '-' nam '-media-memodel.fig']);
saveas(fig(1), ['./../figures/tiff/' model_nam '-' nam '-media-memodel.tif']);

%% Gene Deletion Analyis
unqgenes = unique(model.genes); %all genes in the model
for i = 2:length(unqgenes)    
    modeltemp = model ;
    modeltemp = deleteModelGenes(modeltemp, unqgenes(i));
    model3 = modeltemp;
        
    [solf.x, sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    if ~isempty(solf.x) & ~isnan(solf.x)
        screen_rpmi(i,2) = solf.x(objpos); % impact on growth
    end
     
    model3.c(rxnpos1) = epsilon;
    [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    if ~isempty(solf.x) & ~isnan(solf.x)
        screen_rpmi(i,1) = solf.x(rxnpos1); % impact on flux
    end
    disp(i)
end
    
screen_rpmi(1,:) = NaN;
[sx, spos] = sort(screen_rpmi(:,1));
genessorted = unqgenes(spos);

%% Media culture changes
glucflag = logical([1 0 1 1 1 1 0 0 0 0]);
aaflag =  logical([1 0 1 0 1 0 0 1 0 1]);
pyrflag = logical([1 0 1 1 0 0 0 0 1 1]);
glnflag = logical([1 0 0 1 1 0 1 0 0 1]);
expval = [1 0 1 1 1 1 1 0 1 1];

% part 2
glucflag1 = [ ones(1,6),  0 ,0, 0,  0];
expval1 = [ ones(1,7),  0 , 1,  0];


[ix, aapos] = ismember( mediareactions1(20:37), model.rxns);
posgluc = 1385;  % glucose uptake reaction in recon1.
glnpos = 1386; % glutamine
pyrpos = 1509; % pyruvate
acetatepos = 1238; % acetate
fattyacidpos = 1445; % linoleic acid

model3 = model;
model3.c(rxnpos1) = epsilon;
[solf.x, sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
wild_type = solf.x(rxnpos1); % default flux

for i = 1:10
    model3 = model;
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
    
    model3.c(rxnpos1) = epsilon;
    [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    media_screen_dmem(i,1) = solf.x(rxnpos1); % impact on acetylation flux
    

end
disp('consistency with 10 different experimental conditions - Part I')
sum((media_screen_dmem > 0.05*wild_type) == expval') % 9

for i = 1:10
    model3 = model;
    if ~glucflag1(i)
        model3.lb(posgluc) = 0;
    end
   
    switch i 
        case 2 % no ca+
                    model3.lb(ismember( {'EX_ca2(e)'}, model.rxns)) = 0;
        case 4 % no Phosphate
                    model3.lb(ismember( {'EX_pi(e)'}, model.rxns)) = 0;
        case 6 % no vitamins in dmem.. compmosition from sigma 
                    model3.lb(ismember( {'EX_thm(e)';'EX_ribflv(e)';'EX_pydx(e)';'EX_ncam(e)';'EX_inost(e)';'EX_5mthf(e)';'EX_pnto_R(e)';'EX_chol(e)'},acetylation_model.rxns)) = 0;
        case 7 %  add acetate
                    model3.lb(acetatepos) = -20;
        case 9 % add fatty acid
                model3.lb(fattyacidpos) = -5;
    end

        
    model3.c(rxnpos1) = epsilon;
    [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);
    media_screen_dmem1(i,1) = solf.x(rxnpos1); % impact on flux
end

disp('consistency with 10 different experimental conditions - Part II')
sum((media_screen_dmem1 > 0.05*wild_type) == expval1') % 10

%% Soecific histone markers
for i = 1:14
    % match cell line in CCLE data
    iii = find(ismember(celllinenames_ccle1, acetlevellist(i)));
    if ~isempty(iii)
        iii  = iii(1);
        model2 = model;
        
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
            
        model2.c(rxnpos1) = epsilon;
        [fluxstate_gurobi] =  constrain_flux_regulation(model2, onreactions,...
            offreactions, kappa, rho, epsilon, MODE,[], minfluxflag);
        grate_ccle_exp_acetdat(i,2) = fluxstate_gurobi(rxnpos1); %acetylation flux    
    end
 end
 
figure;
plot(grate_ccle_exp_acetdat(:,2), acet_meth_listval(8, :)', 'o',...
    'markerfacecolor', [0.9020, 0.3804, 0.0039], 'markeredgecolor', 'k') 
fb = polyfit(grate_ccle_exp_acetdat(:,2), acet_meth_listval(8,:)',1);
fb1 = polyval(fb, grate_ccle_exp_acetdat(:,2));

hold on; 
plot(grate_ccle_exp_acetdat(:,2), fb1, 'r-', 'linewidth', 2);
text(grate_ccle_exp_acetdat(:,2) + 0.2, acet_meth_listval(8, :)', acetlevellist,...
    'Margin', 5);%, 'editing', 'on');
    %'horizontalalignment', 'right', 'fontsize', 10, 'fontname', 'helvetica')
    %'fontweight', 'bold')
xlabel('SAM flux (mmol/gDW*hr)');
ylabel('H3K9me1');

[acetlevelcorr, acetlevelcorrpv] = corr(grate_ccle_exp_acetdat(:,2), acet_meth_listval(8,:)');


