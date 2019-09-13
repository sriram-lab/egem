function group_pval = fisher_pvalue_meta_analysis(pvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute Fisher's p-value for meta-analysis
%    pvals is a vector which holds p-values
%
%    Fisher's method of combining p-values described in
%       https://en.wikipedia.org/wiki/Fisher%27s_method
%
%    Code adapted from 
%       https://stats.stackexchange.com/questions/158225/why-is-my-combined-p-value-obtained-using-the-fishers-method-so-low
%       credit user Dmitry Smirnov
%           https://stats.stackexchange.com/users/75214/dmitry-smirnov
%
%    Example usage:
%       final_p_value = fisher_pvalue_meta_analysis([0.1,0.2])
%       final_p_value = fisher_pvalue_meta_analysis([0.01,0.02,0.003])
%       final_p_value = fisher_pvalue_meta_analysis([0.01,0.02,0.003, 0.001, 0.01, 0.04])
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pvals is a vector which holds p-values
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals))
    nsig = sum(pvals < 0.05);
    
    hist(pvals)
    title("Histogram of p-values")
    xlabel("p-values"); ylabel("Frequency")
    
    