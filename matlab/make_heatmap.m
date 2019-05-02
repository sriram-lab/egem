%% make_heatmap.m
function data = make_heatmap(model, compartment, epsilon2)
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
load conditions, metabolites

% Initial model parameters
minfluxflag = 0;
posgluc = 1385;  % glucose uptake reaction in RECON1
objpos = find(model.c); % biomass objective
epsilon2 = 1; % higher weights for methylation compared to acetylation

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
        
        % add demand reactions from the metabolite list
        for m=1:length(metabolites)
            tmp_met = char(metabolites(m,2));
            tmp = [tmp_met,'[',compartment,'] -> '];
            tmpname = char(metabolites(m,1));
            model2 = addReaction(model2, tmpname,...
                'reactionFormula', tmp);
            rxnpos(m,1)  = find(ismember(model2.rxns, metabolites(m,1)));
            name(m,1) = metabolites(m,1);
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
        for rxn=1:length(rxnpos)
            model3.c(rxn) = epsilon2;
            [solf.x, sol11] =  constrain_flux_regulation(model3,[],[],0,0,0,[],[],minfluxflag);

            str = ['media_xchange_',num2str(kappatype),'(j,rxn) = solf.x(rxn);'];
            if ~isempty(solf.x) &  ~isnan(solf.x)
                eval(str)
            end
            disp(j)
        end
        disp(kappatype)
    end
end

% Prepare figure labels and variables
medium_labels = conditions(:,1);
reaction_labels = new(:,3)
excess = media_xchange_1;
depletion = media_xchange_2;

% Heatmap
fig = figure;

subplot(1,2,1);
heatmap(excess)
ax1 = gca;
ax1.XData = reaction_labels;
ax1.YData = medium_labels;

subplot(1,2,2);
heatmap(depletion)
ax2 = gca;
ax2.XData = reaction_labels;
ax2.YData = medium_labels;

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