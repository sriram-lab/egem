MATLAB directory structure:
`/figures`: output from scripts.
  `acetyation`: figures from the acetylation model.
  `Difference`: contains figures generated from taking the difference from a model output with two objective functions versus a model output with one objective function.

`/models`: genome scale metabolic models {MATLAB structure}
  `acety2.mat`: new acetylation model with two additional reactions that convert pyruvate and acetate into acetyl-CoA. This model reduces noise from the original acetylation model.
  `METATn.mat`: a very simplistic methylation model that only includes the conversion of L-methionine + ATP -> SAM via MATII. Does not contains other relevant methylation reactions.
  `model.mat`: a more complex methylation model that includes other methylation reactions.
  `recon1.mat`: the base model. should compare results to this model.

`/new_var`: MATLAB data variables that are not in Shen et al., 2019
  `mem.json`: Minimal medium components + additional ones of interest. Same as Shen et al., 2019.
  `metabolite.json`: map of metabolite identiers
  `reaction_name.json`: map of reaction identifiers.

`/GCP`: Histone chemoproteomics data from the Broad Institute / MIT
  * May be better to use this instead of chemogenomics datasets - need to validate with comparable cell lines from the CCLE

`/scripts`: MATLAB scripts for running analyses
  * `cellline_heatmap`: Get cell-specific responses to metabolic perturbations for specific reactions
  * `histone_corr`: Calculate correlation between histone marker proteomics and metabolic flux
  * `make_heatmap`: Simple heatmap with various demand reaction values in response to media perturbations
  * `read_json`: read json files and make into matlab variables
  * `run_all`: run modules
  * `leroy_names`: the name of the leroy et al 2013 cell lines
  * `metabolites`

`/tables`: Data table outputs from scripts
`/test`: Test directory / sanity checks
