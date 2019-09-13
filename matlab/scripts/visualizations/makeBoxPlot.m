function makeBoxPlot(flux, proteomics, type, folderName, exp)
    addpath('./../visualizations/BrewerMap');
    addpath('./../visualizations/export_fig');
    groups = {'KAC'; 'KMe1'; 'KMe2'; 'KMe3'};
    for rxn = 1:length(groups)
        [HighProteomics, LowProteomics] = stratify(flux(:, rxn), ...
            proteomics, 'mean');

        if isempty(HighProteomics)
            HighProteomics = zeros([5, 1]);
        end
        if isempty(LowProteomics)
            LowProteomics = ones([5, 1]);
        end
        vars = [zeros(size(HighProteomics)); ones(size(LowProteomics))];
        
        subplot(2, 2, rxn)
        boxplot([HighProteomics; LowProteomics], vars, ...
            'orientation', 'vertical', ...
            'symbol', '', ...
            'labels', {'High Flux', 'Low Flux'});
        map = brewermap(2, 'Set1');
        
        boxes = findobj(gca, 'Tag', 'Box');
        
        for boxNum = 1:length(boxes)
            jj = length(boxes) - boxNum + 1;
            patch(get(boxes(boxNum), 'XData'), ...
                get(boxes(boxNum), 'YData'), ...
                map(jj,:), ...
                'FaceAlpha', 0.20, ...
                'edgecolor', 'k');
        end

        set(gca, 'box', 'off')
        set(gca, 'linewidth', 1.3)
        set(gcf, 'color', 'white')
        grid('on')

        set(gca, 'fontsize', 11)
        set(findobj(gca, 'Type', 'text'), 'FontSize', 11)
        set(findobj(gca, 'Type', 'line'), 'linewidth', 1.2)
        set(gca,'TickDir', 'out')
        ylabel('PTM Levels')
        title(string(groups(rxn)), 'FontSize', 11)
        set(gcf, 'position', [470.0000, 589.5000, 341.5000, 308.5000])

        base = './../../../figures/';
        strs = strcat(base, folderName, '/', type, exp, 'BP');
        %export_fig(strs, '-r250', '-tif');
        strs1 = strcat(strs, '.fig');
        saveas(gcf, strs1);
    end
end