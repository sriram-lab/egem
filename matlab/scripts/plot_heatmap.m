%% plot_heatmap.m 
% Author: Scott Campit
function plot_heatmap(STRUCT,...
    metabolites, exp, epsilon2, medium, type)

%% Generate Heatmaps for Downstream Analysis

load('./../vars/mediareactions.mat') % Medium components
load('./../vars/metabolites.mat') % Reactions you're interested in observing

switch exp
    case {'sra', 'competition', 'no_competition'} 
        if ~isequal(size(epsilon2), [1,1])
            epsilon2 = '';
        end
        medium_labels = mediareactions(:,2);
        reaction_labels = metabolites(:,3);

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_flux_rc))
        ax6 = gca;
        ax6.XData = reaction_labels;
        ax6.YData = medium_labels;
        ax6.Title = 'Reduced cost in depleted medium';
        xlabel(ax6, 'Demand reactions');
        ylabel(ax6, 'Medium component');
        
        % Growth rate
        fig4 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_grate))
        ax5 = gca;
        ax5.XData = reaction_labels;
        ax5.YData = medium_labels;
        ax5.Title = 'Growth rate in excess medium';
        xlabel(ax5, 'Demand reactions');
        ylabel(ax5, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_grate))
        ax6 = gca;
        ax6.XData = reaction_labels;
        ax6.YData = medium_labels;
        ax6.Title = 'Growth rate in depleted medium';
        xlabel(ax6, 'Demand reactions');
        ylabel(ax6, 'Medium component');
        
        base = strcat('./../figures/new-model/eGEMn_',...
            string(epsilon2),...
            '_', exp, medium);
        
        % Create filenames for all figures
        fig1_str = strcat(base, '_flux.fig');
        fig2_str = strcat(base, '_sp.fig');
        fig3_str = strcat(base, '_rc.fig');
        fig4_str = strcat(base, '_grate.fig');

        saveas(fig1, char(fig1_str));
        saveas(fig2, char(fig2_str));
        saveas(fig3, char(fig3_str));
        saveas(fig4, char(fig4_str));
    
        base = strcat('./../figures/new-model/', ...
            medium, '_', ...
            string(epsilon2), '_', ...
            exp);
        
        % Create filenames for all figures
        fig1_str = strcat(base, '_flux.fig');
        fig2_str = strcat(base, '_shadowPrice.fig');
        fig3_str = strcat(base, '_reducedCosts.fig');
        %fig4_str = strcat(base, '_grate.fig');

        saveas(fig1, fig1_str);
        saveas(fig2, fig2_str);
        saveas(fig3, fig3_str);
        %saveas(fig4, char(fig4_str));
    
    case 'fva'  
        medium_labels = mediareactions(:,2);
        reaction_labels = metabolites(:,3);

        fig = figure;
        subplot(1,2,1);
        heatmap(normalize(STRUCT.excess_maxflux))
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(STRUCT.depletion_maxflux))
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        base = strcat('./../figures/new-model/',...
            medium, '_', ...
            epsilon2, '_', ...
            exp);
        
        % Create filenames for all figures
        fig_str = strcat(base, '.fig');
        saveas(fig, char(fig_str));
    
    case 'correlation'
        
        % Get data
        rho = STRUCT.R;
        rxns = STRUCT.Reaction;
        marks = STRUCT.HistoneMark;
        
        fig = figure;
        
        heatmap(rho)
        ax = gca;
        ax.XData = marks;
        ax.YData = rxns;
        ax.Title = 'Histone markers and metabolic flux correlation';
        xlabel(ax, 'Histone Markers');
        ylabel(ax, 'Cancer Cell Lines');
        str = strcat('./../figures/corr/', STRUCT.Name, type, '.fig');
        saveas(fig, str);
    
    case 'pval'
        % Get data and set values that are not significant to NaN
        pval = STRUCT.Pvalue;
        rmv = (pval > 0.05);
        pval(rmv) = 0;
        pval(pval == 0) = NaN;
        
        rxns = STRUCT.Reaction;
        marks = STRUCT.HistoneMark;
        
        fig = figure;
        
        heatmap(pval)
        ax = gca;
        ax.XData = marks;
        ax.YData = rxns;
        ax.Title = 'Histone markers and metabolic flux correlation';
        xlabel(ax, 'Histone Markers');
        ylabel(ax, 'Demand Reactions');
        str = strcat('./../figures/corr/', STRUCT.Name, type, '.fig');
        saveas(fig, str);
        % Then use Plotly offline
%         fig = plotlyfig(gcf);
%         fig.layout.showlegend = true;
%         fig.layout.title = 'Correlation between LeRoy histone markers and eGEM flux';
% 
%         fig.PlotOptions.FileName = 'leroy_histone_corr';
%         fig.PlotOptions.FielOpt = 'overwrite';
% 
%         getplotlyoffline('https://cdn.plot.ly/plotly-latest.min.js') 
%         fig.PlotOptions.Offline = true;
% 
%         fig.plotly
end
end