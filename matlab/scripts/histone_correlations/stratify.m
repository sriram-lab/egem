function [HighProteomics, LowProteomics] = stratify(flux, proteomics, thresh)
    proteomics = normalize(proteomics, 'range');
    switch thresh
        case 'median'
            threshold = median(flux, 'all');
            HighProteomics = proteomics(flux > threshold);
            LowProteomics = proteomics(flux < threshold);
        case 'mean'
            threshold = mean(flux, 'all');
            HighProteomics = proteomics(flux > threshold);
            LowProteomics = proteomics(flux < threshold);
    end
end