function plotHistoneCorrelation(STRUCT, var, type, exp)
    switch var
        case 'correlation'
            % Get data
            rho = STRUCT.PearsonR;
            rxns = STRUCT.Reaction;
            marks = STRUCT.HistoneMark;

            fig = figure;

            heatmap(rho)
            ax = gca;
            ax.XData = marks;
            ax.YData = rxns;
            xlabel(ax, 'Histone Markers');
            ylabel(ax, 'Reactions');
            
            switch exp
                case 'tissue_corr'
                    ax.Title = string(strcat(string(STRUCT.Tissue),...
                        ' histone markers and metabolic flux correlation'));
                    str = strcat('./../../../figures/tissue_corr/', ...
                        STRUCT.Tissue, type, 'Corr.fig');
                case 'medium_corr'
                    ax.Title = string(strcat(string(STRUCT.Medium),...
                        ' histone markers and metabolic flux correlation'));
                    str = strcat('./../../../figures/medium_corr/', ...
                        STRUCT.Medium, type, 'Corr.fig');
                case 'culture_corr'
                    ax.Title = string(strcat(string(STRUCT.Culture),...
                        ' histone markers and metabolic flux correlation'));
                    str = strcat('./../../../figures/culture_corr/', ...
                        STRUCT.Culture, type, 'Corr.fig');
            end
            saveas(fig, str);

        case 'pval'
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
            ylabel(ax, 'Reactions');
            str = strcat('./../../../figures/culture_corr/', STRUCT.Culture, 'Corr.fig');
            saveas(fig, str);
    end
end