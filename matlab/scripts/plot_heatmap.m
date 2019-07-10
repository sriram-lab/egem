%% plot_heatmap.m 
% Author: Scott Campit
function plot_heatmap(struct, exp, epsilon2, medium)

%% Generate Heatmaps for Downstream Analysis
load ./../vars/supplementary_software_code media_exchange1
var = {'./../vars/metabolites.mat', './../vars/cellmedia.mat',...
    './../vars/mediareactions1.mat'};
for kk = 1:numel(var)
    load(var{kk});
end

switch exp
    case 'sra' 
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);
        
        % Metabolic fluxes
        fig1 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux))
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_flux))
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        % Shadow prices
        fig2 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux_sp))
        ax3 = gca;
        ax3.XData = reaction_labels;
        ax3.YData = medium_labels;
        ax3.Title = 'Shadow price in excess medium';
        xlabel(ax3, 'Demand reactions');
        ylabel(ax3, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_flux_sp))
        ax4 = gca;
        ax4.XData = reaction_labels;
        ax4.YData = medium_labels;
        ax4.Title = 'Shadow price in depleted medium';
        xlabel(ax4, 'Demand reactions');
        ylabel(ax4, 'Medium component');
        
        % Reduced cost
        fig3 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux_rc))
        ax5 = gca;
        ax5.XData = reaction_labels;
        ax5.YData = medium_labels;
        ax5.Title = 'Reduced cost in excess medium';
        xlabel(ax5, 'Demand reactions');
        ylabel(ax5, 'Medium component');

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

        saveas(fig1, fig1_str);
        saveas(fig2, fig2_str);
        saveas(fig3, fig3_str);
        saveas(fig4, fig4_str);
    
    case 'fba' 
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);
        
        % Metabolic fluxes
        fig1 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux))
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_flux))
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        % Shadow prices
        fig2 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux_sp))
        ax3 = gca;
        ax3.XData = reaction_labels;
        ax3.YData = medium_labels;
        ax3.Title = 'Shadow price in excess medium';
        xlabel(ax3, 'Demand reactions');
        ylabel(ax3, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_flux_sp))
        ax4 = gca;
        ax4.XData = reaction_labels;
        ax4.YData = medium_labels;
        ax4.Title = 'Shadow price in depleted medium';
        xlabel(ax4, 'Demand reactions');
        ylabel(ax4, 'Medium component');
        
        % Reduced cost
        fig3 = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_flux_rc))
        ax5 = gca;
        ax5.XData = reaction_labels;
        ax5.YData = medium_labels;
        ax5.Title = 'Reduced cost in excess medium';
        xlabel(ax5, 'Demand reactions');
        ylabel(ax5, 'Medium component');

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
            '_', exp, '_', medium);
        
        % Create filenames for all figures
        fig1_str = strcat(base, '_flux.fig');
        fig2_str = strcat(base, '_sp.fig');
        fig3_str = strcat(base, '_rc.fig');
        fig4_str = strcat(base, '_grate.fig');

        saveas(fig1, fig1_str);
        saveas(fig2, fig2_str);
        saveas(fig3, fig3_str);
        saveas(fig4, fig4_str);
    
    case 'fva'  
        medium_labels = mediareactions1(:,2);
        reaction_labels = metabolites(:,3);

        fig = figure;
        subplot(1,2,1);
        heatmap(normalize(struct.excess_max_flux))
        ax1 = gca;
        ax1.XData = reaction_labels;
        ax1.YData = medium_labels;
        ax1.Title = 'Metabolic flux in excess medium';
        xlabel(ax1, 'Demand reactions');
        ylabel(ax1, 'Medium component');

        subplot(1,2,2);
        heatmap(normalize(struct.depletion_max_flux))
        ax2 = gca;
        ax2.XData = reaction_labels;
        ax2.YData = medium_labels;
        ax2.Title = 'Metabolic flux in depleted medium';
        xlabel(ax2, 'Demand reactions');
        ylabel(ax2, 'Medium component');
        
        base = strcat('./../figures/new-model/eGEMn_', exp, '_', medium);
        
        % Create filenames for all figures
        fig_str = strcat(base, '_flux.fig');
        saveas(fig, fig_str);

end
end