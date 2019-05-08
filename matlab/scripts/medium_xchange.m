%% medium_xchange
function medium_xchange()
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
                rxn = rxn;
            else
                model2 = addReaction(model2, ['DM_', metabolite], 'reactionFormula', [metabolite '[' compartment '] -> ']);
                rxn  = [find(ismember(model2.rxns, ['DM_' metabolite]));];
                nam = ['DM_' metabolite];
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
            model3.c(rxn) = new_epsilon;
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