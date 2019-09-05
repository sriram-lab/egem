function makeHist(flux, name, type, exp)
    
    bins = 100;

    fig = figure;
    h1 = histogram(flux(:, 1), bins);
    hold on;
    h2 = histogram(flux(:, 2),bins);
    hold on;
    h3 = histogram(flux(:, 3), bins);
    hold on;
    h4 = histogram(flux(:, 4), bins);
    set(gca, 'yscale', 'log');
    
    title("Distribution of flux among CCLE cell lines");
    xlabel("Predticted flux");
    ylabel("Total number of cell lines");
    legend("Ac", "1Me", "2Me", "3Me");
    
    str = strcat('./../../../figures/', type, '/', name, exp, 'Hist.fig');
    saveas(fig, str);
end