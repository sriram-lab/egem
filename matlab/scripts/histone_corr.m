%% histone_corr calculates the correlation value between various histone markers and the metabolic flux obtained from the iMAT algorithm
function histone_corr(model, metabolite, rxn, compartment, mode, epsilon, rho, kappa, minfluxflag)
%% INPUTS:
    % model: Initial genome scale model
    % metabolite: The metabolite we're interested in creating a demand
    % reaction. If this argument is left empty, it will take the rxn
    % argument.
    % rxn: The rxn of interest that we are creating a demand reaction through.
    % compartment: specifying the compartment. 
    % mode: for constrain flux regulation
    % epsilon: for constrain flux regulation
    % rho: 
    % kappa:
    % minfluxflag:
%% OUTPUTS:
    % histogram plotting the correlation value (x-axis) corresponding
    % to each histone marker (y-axis) for the demand reaction of
    % interest.
%% histone_corr
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz %contains CCLE cell line names, gene expression data (z-transformed)
load supplementary_software_code acetlevlistmedia acetlevellist acetlevellistval %contains cell line names, growth media , total bulk acetylation

path = './../new_var/';
vars = {[path 'h3marks.mat'], [path 'h3names.mat'], [path 'h3vals.mat']};



for kk = 1:numel(var)
    load(var{kk})
end

iii = find(ismember(celllinenames_ccle1, ));
if ~isempty(iii)
    iii  = iii(1);
    model2 = model;

    % Takes in genes that are differentially expression from Z-score
    % scale
    ongenes = unique(ccleids_met(ccle_expression_metz(:,iii) >= 2));
    offgenes = unique(ccleids_met(ccle_expression_metz(:,iii) <= -2));

    % now set the media and glucose levels for different media conditions
    if ismember({'RPMI'} , acetlevlistmedia(i)) % RPMI
        model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5;
    elseif ismember({'DMEM'} , acetlevlistmedia(i)) % DMEM
        model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*4.5/2;
    elseif ismember({'L15'} , acetlevlistmedia(i)) % NO GLUC AND LOW GAL  
        model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -0; % L15 
        model2.lb(find(ismember(model2.rxns, {'EX_gal(e)'}))) = -0.9; % LOW GAL
    elseif ismember({'McCoy 5A'} , acetlevlistmedia(i)) % McCoy
        model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*3/2;
    elseif ismember({'IMM'} , acetlevlistmedia(i)) % IMDM
        model2.lb(find(ismember(model2.rxns, {'EX_glc(e)'}))) = -5*4.5/2;
    end

    % Get the reactions corresponding to on- and off-genes
    [~,~,onreactions,~] =  deleteModelGenes(model2, ongenes);
    [~,~,offreactions,~] =  deleteModelGenes(model2, offgenes);
    disp(i)

    % Get the flux redistribution values associated with different media component addition and deletion
    [fluxstate_gurobi, grate_ccle_exp_dat(i,1),  solverobj_ccle(i,1)] =  constrain_flux_regulation(model2,onreactions,offreactions,kappa,rho,epsilon,MODE ,[], minfluxflag);

    % Now let's add the demand reaction we want for a given metabolite
    if (~exist('metabolite', 'var')) || (isempty(metabolite))
        rxn = rxn;
    else
        model2 = addReaction(model2, ['DM_', metabolite], 'reactionFormula', [metabolite '[' compartment '] -> ']);
        rxn  = [find(ismember(model2.rxns, ['DM_' metabolite]));];
        nam = ['DM_' metabolite];
    end

    % limit methionine levels for all reactions in the model; it has to be non limiting
    [ix, pos]  = ismember({'EX_met_L(e)'}, model2.rxns);
    model2.lb(pos) = -0.5;
    model2.c(3743) = 0;
    model2.c(rxn) = 0.01; % we're interested in this reaction

    % get the flux values from iMAT
    [fluxstate_gurobi] =  constrain_flux_regulation(model2,...
        onreactions, offreactions, kappa, rho, epsilon, MODE ,[],...
        minfluxflag);
    grate_ccle_exp_dat(i,2) = fluxstate_gurobi(rxn);
end

% Calculate the pearson correlation coefficient
[acetlevelcorr_amet, acetlevelcorrpv_amet] = corr(grate_ccle_exp_dat(:,2),...
    acet_meth_listval');
acetlevelcorr_amet = acetlevelcorr_amet';

% Make plot of correlation coefficients versus histone markers
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
saveas(fig(1), ['./../figures/fig/methylation-' nam '-corr.fig']);
saveas(fig(1), ['./../figures/tiff/methylation-' nam '-corr.tif']);
end 