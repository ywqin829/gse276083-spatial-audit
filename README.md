# Spatial attribution of candidate expression gradients in GSE276083

This exploratory reanalysis asks whether gradients around cells with neural-like transcriptomic labels can also be recovered around nearby non-neural cells. Harutaka Katano and Yoshihiko Hoshino contributed the original [GEO series GSE276083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276083). The archived GEO count, annotation and coordinate tables are the source team's data. The section selection, revised labels, distances, model fits and reference comparisons are products of the secondary analysis.

## Reproduce the audit

Extract analysis_code.zip into a folder named code, next to section_inputs.rds. From code, run the commands below. The primary run takes the November reference seed 42. The sensitivity run uses seeds 1 through 20. Tested with R 4.5.1, data.table 1.17.8, Matrix 1.7-3 and RANN 2.6.2.

```sh
Rscript --vanilla run_audit.R ../section_inputs.rds candidates.csv ../reproduced 42
Rscript --vanilla seed_sensitivity.R ../section_inputs.rds candidates.csv ../reproduced
```

The RDS holds ordered cell identifiers, section labels, broad and fine cell annotations, global coordinates in pixels, counts, segmented area, negative-probe counts per probe, neighbourhood covariates, cached reference distance and normalized expression for 19 unique genes. The candidate list specifies ten gene/cell-type pairs per disease and includes CXCL16 in both groups. The starting point for exact reproduction is this prepared input. The included export_inputs.R documents extraction from the November reannotated Seurat object, distance cache and neighbourhood table using Seurat 5.3.0. Reconstructing the original manual annotations from the GEO tables alone is outside this release.

## Sections and models

The six-section analysis includes 337,875 cells from three MTB and three NTM lung sections. Section 18_77_Bronch, identified as a bronchial specimen in the secondary analysis, was excluded by section ID. The RDS preserves the cell-level inputs while all summaries give equal weight to sections. The seven-row inclusion table records the original group field separately from the analysis classification. NTM species were unavailable for species-specific comparisons.

In each section and target cell type, expression is regressed on distance from the nearest neural-like cell. Slopes are per 100 pixels. M0 includes distance and an intercept. M3 adds log-transformed total counts, cell area and neighbour count, negative-probe counts per probe and four neighbourhood composition covariates. The latter are Epithelial, Immune, Stromal and Unclassified. Models use the eligible complete cases for their design and omit invariant adjustment columns. Retention is the mean M3 slope divided by the mean M0 slope across three sections. Ratios are unstable for near-zero baseline slopes. The retained-slope rule requires a positive adjusted slope in three of three sections and retention of at least 0.40.

For each neural-like anchor the position control samples one non-neural cell from up to 25 nearest neighbours within 100 pixels. Matches may be reused across anchors. Unmatched anchors contribute no match and duplicate references are collapsed before distance calculation. A mast-cell control uses cells labelled Mast_cells in the same section. The distance calculation excludes self-matches by cell index. The code preserves the random-number sequence of 20 abundance-matched draws between sections used in the November analysis. Section-wise Spearman correlations describe expression versus distance for each reference. A reproduced candidate direction means positive correlations in all three sections. The matching sensitivity analysis holds the tissue and expression fixed and varies only the selected non-neural cells. These draws give no additional biological replicates.

## Reading the results

| File in audit_results.zip | Contents |
| --- | --- |
| S0_section_cell_counts.csv | Annotated cell counts and inclusion flag per section |
| section_estimates.csv | M0 and M3 slopes, fitted cell counts and three correlations per candidate and section |
| candidate_summary.csv | Section counts, retention and reference contrasts for seed 42 |
| reference_summary.csv | Number of candidates with positive correlations in three of three sections, by reference and disease |
| reference_qc.csv | Reference sizes, maximum matching offsets and comparison with the November distance cache |
| matched_pairs.csv | Anchor IDs, selected non-neural IDs and offsets for the primary draw |
| seed_candidates.csv | Candidate-level statistics for each of 20 sensitivity draws |
| seed_summary.csv | Counts by disease and seed |

At seed 42, position-matched references reproduced 8 of 10 MTB and 9 of 10 NTM candidate directions. The respective mast-cell counts were 1 and 0. Across seeds 1 to 20, the position-matched counts ranged from 8 to 10 in each group. In each NTM section epithelial ANXA1, GPX3 and MT2A had positive M3 slopes, with group-level retention of 0.759, 0.699 and 0.511, respectively.

The field specific is an operational threshold on the correlation contrast. It requires a positive neural-minus-matched difference in all three sections and a mean difference above 0.01. The stricter sensitivity field specific_all_above_margin requires a difference above 0.01 in each section. At seed 42, NTM ITGA5 and MT2A satisfy both the retained-slope rule and the mean-margin rule. Only MT2A also satisfies the stricter per-section margin. Across the 20 sensitivity draws no candidate consistently satisfies the combined retained-slope and mean-margin rule. These thresholds describe reference sensitivity. Direct evidence for a nerve-mediated effect would require independent anatomical identification and functional evidence. The selected candidates and three-section groups limit generalization.
