function [stats] = DescriptiveStats(array)
    for col = 1:size(array, 2)
        allMax(1, col) = max(array(:, col), 'omitnan');
        allMean(1, col) = mean(array(:, col), 'omitnan');
        allMedian(1, col) = median(array(:, col), 'omitnan');
        allMin(1, col) = min(array(:, col), 'omitnan');
        allStd(1, col) = std(array(:, col), 'omitnan');
        allVar(1, col) = var(array(:, col), 'omitnan');
    end
    stats.Max = allMax;
    stats.Mean = allMean;
    stats.Median = allMedian;
    stats.Min = allMin;
    stats.Std = allStd;
    stats.Var = allVar;
end