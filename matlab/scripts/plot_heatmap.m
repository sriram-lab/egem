%% plot_heatmap.m 
% Author: Scott Campit
function plot_heatmap(STRUCT,...
    metabolites, exp, epsilon2, medium, type)

%% Generate Heatmaps for Downstream Analysis
<<<<<<< HEAD
load('./../vars/mediareactions.mat') % Medium components
load('./../vars/metabolites.mat') % Reactions you're interested in observing
=======
load ./../vars/supplementary_software_code media_exchange1
var = {...
    './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'...
    };
for kk = 1:numel(var)
    load(var{kk});
end
>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967

switch exp
    case {'sra', 'competition', 'no_competition'} 
        if ~isequal(size(epsilon2), [1,1])
            epsilon2 = '';
        end
<<<<<<< HEAD
        medium_labels = mediareactions(:,2);
        reaction_labels = metabolites(:,3);

=======
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);
        
>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967
        % Metabolic fluxes
        fig1 = figure;
        subplot(1,2,1);
        heatmap(normalize(STRUCT.excess_flux))
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(STRUCT.depletion_flux))
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        % Shadow prices
        fig2 = figure;
        subplot(1,2,1);
        heatmap(normalize(STRUCT.excess_flux_sp))
        ax3 = gca;
        ax3.XData = reaction_labels;
        ax3.YData = medium_labels;
        ax3.Title = 'Shadow price in excess medium';
        xlabel(ax3, 'Demand reactions');
        ylabel(ax3, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(STRUCT.depletion_flux_sp))
        ax4 = gca;
        ax4.XData = reaction_labels;
        ax4.YData = medium_labels;
        ax4.Title = 'Shadow price in depleted medium';
        xlabel(ax4, 'Demand reactions');
        ylabel(ax4, 'Medium component');
        
        % Reduced cost
        fig3 = figure;
        subplot(1,2,1);
        heatmap(normalize(STRUCT.excess_flux_rc))
        ax5 = gca;
        ax5.XData = reaction_labels;
        ax5.YData = medium_labels;
        ax5.Title = 'Reduced cost in excess medium';
        xlabel(ax5, 'Demand reactions');
        ylabel(ax5, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(STRUCT.depletion_flux_rc))
        ax6 = gca;
        ax6.XData = reaction_labels;
        ax6.YData = medium_labels;
        ax6.Title = 'Reduced cost in depleted medium';
        xlabel(ax6, 'Demand reactions');
        ylabel(ax6, 'Medium component');
        
%         % Growth rate
%         fig4 = figure;
%         subplot(1,2,1);
%         heatmap(normalize(struct.excess_grate))
%         ax5 = gca;
%         ax5.XData = reaction_labels;
%         ax5.YData = medium_labels;
%         ax5.Title = 'Growth rate in excess medium';
%         xlabel(ax5, 'Demand reactions');
%         ylabel(ax5, 'Medium component');
% 
%         subplot(1,2,2);
%         heatmap(normalize(struct.depletion_grate))
%         ax6 = gca;
%         ax6.XData = reaction_labels;
%         ax6.YData = medium_labels;
%         ax6.Title = 'Growth rate in depleted medium';
%         xlabel(ax6, 'Demand reactions');
%         ylabel(ax6, 'Medium component');
        
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
<<<<<<< HEAD
        medium_labels = mediareactions(:,2);
=======
        medium_labels = mediareactions1(:,2);
>>>>>>> c9c95c7ea86416cfdce0cfebfda8aab4fd025967
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