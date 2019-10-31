function checkMeanCentering(array)
    for predictor = 1:size(array, 2)
        hist(array(:, predictor), 1000)
    end
end

