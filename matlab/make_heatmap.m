%% make_heatmap.m
function data = make_heatmap(model, rxnpos, epsilon2)
% make_heatmap.m displays an array of correlation values corresponding to
% several demand reactions and the quantity of the histone maker based on
% metabolic perturbation. 

% INPUTS:
    % model: genome scale metabolic model
    % rxnpos: list of reactions that will be used in the analysis
    % epsilon2: a weight that will be used for the secondary objective coefficient 
    
% load datasets into function workspace
load supplementary_software_code acetylation_model 
load supplementary_software_code celllinenames_ccle1 ccleids_met ccle_expression_metz 
load supplementary_software_code cellmedia acetlevellist acetlevellistval 
load methylation_proteomics_validation_data acet_meth_listval acet_meth_list_rowlab 
load supplementary_software_code labels media_exchange1 mediareactions1
load conditions

% Initial model parameters
minfluxflag = 0;
posgluc = 1385;  % glucose uptake reaction in RECON1
objpos = find(model.c); % biomass objective
epsilon2 = 1; % higher weights for methylation compared to acetylation

%% First calculate t
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
    
% Kappa values are weights for both excess and depletion of medium.
% kappa=10 is excess, and kappa=0.01 is depletion
for kappatype = 1:2
    if kappatype == 1 
        kappa  = 10; 
    else
        kappa = 0.01;
    end

    % Different kappa values will be assigned for trace minerals and
    % biomass components relative to the other medium components.
    for j = 1:50
        kappa1 = kappa;
        if (kappatype == 2) & (ismember(j,[1:7, 28:35, 37, 38, 43])) % trace elements: lower kappa
            kappa1 = kappa/100;
        elseif (kappatype == 1) & (ismember(j,[2, 13, 36])) % glucose or glutamine: higher kappa
            kappa1 = 3;
        end
        model2 = model;

        % Make the methylation exchange reaction have a fixed LB of
        % -0.5 to be non-limiting
        [ix, pos]  = ismember({'EX_met_L(e)'}, model2.rxns);
        model2.lb(pos) = -0.5; 

        % Identify the reaction of interest. 
        if (~exist('meth_type','var')) || (isempty(meth_type))
            rxnpos = rxnpos;
        else
            model2 = addReaction(model2, ['DM_', metabolite], 'reactionFormula', [metabolite '[' compartment '] -> ']);
            rxn  = [find(ismember(model2.rxns, ['DM_' metabolite]));];
            nam = ['DM_' metabolite];
        end

        % Get the reaction position for the various medium components
        % and constrain the lower bound with a value proportional to
        % parameter kappa1
        [ix, pos]  = ismember(conditions(j,2), model2.rxns);
        model2.lb(pos) = -1*cell2mat(conditions(j,3))*kappa1;

        % Obtain the flux values in excess media conditions using the iMAT 
        % algorithm with reaction lower bound values that are proportional
        % to kappa1.
        [solf.x, sol11] =  constrain_flux_regulation(model2,[],[],0,0,0,[],[],minfluxflag);

        str = ['media_change_growth_', num2str(kappatype), '(j,1) = solf.x(objpos);'];
        if ~isempty(solf.x) & ~isnan(solf.x)
            eval(str)
        end

        % Obtain flux values when using epsilon2 as the objective
        % coefficient for the reaction of interest.
        model3 = model2;
        model3.c(rxnpos) = epsilon2;
        [solf.x,sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);

        str = ['media_xchange_',num2str(kappatype),'(j,1) = solf.x(rxnpos);'];
        if ~isempty(solf.x) &  ~isnan(solf.x)
            eval(str)
        end

        disp(j)
    end
    disp(kappatype)
end

for x = 1:length(labels)
    [R_up, pval_up] = corr(ccle_grate(x,2), media_xchange_1(:, 1));
    [R_down, pval_down] = corr(ccle_grate(x,2), media_xchange_2(:, 1));

% Make the figure
labels = conditions(:,1);
labels(2) = [];

fig = figure;
data = [media_xchange_1(:, 1), media_xchange_2(:, 1)];
data(2,:) = [];
plt = barh(data, 'edgecolor', 'w');
set(plt(2), 'FaceColor', hex2rgb('#C6393D'));
set(plt(1), 'FaceColor', hex2rgb('#BDCD6C'));
%title('Varying media components', 'fontweight', 'bold');
set(gca,'ytick', [1:length(labels)], 'yticklabel',...
    labels, 'fontsize', 8, 'fontweight', 'bold');
set(gca,'TickDir', 'out');
set(gca,'box', 'off');
set(gca,'linewidth', 2);
set(gcf,'color', 'white');
set(gca,'fontsize', 12);
set(gcf, 'Position', [100, 100, 700, 800])
xlabel('Demand reaction flux (mmol/gDW*hr)');
ylabel('Medium components')
h = legend({'Excess', 'Depletion'});
legend boxoff;
%saveas(fig(1), ['./../figures/fig/' model_nam '-' nam '-media-memodel.fig']);
%saveas(fig(1), ['./../figures/tiff/' model_nam '-' nam '-media-memodel.tif']);
end