# Spatial attribution of candidate expression gradients in GSE276083

This exploratory reanalysis asks whether gradients around cells with neural-like transcriptomic labels can also be recovered around nearby non-neural cells. Harutaka Katano and Yoshihiko Hoshino contributed the original [GEO series GSE276083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276083). The archived GEO count, annotation and coordinate tables are the source team's data. The section selection, revised labels, distances, model fits and reference comparisons are products of the secondary analysis.

The version-controlled analysis code is maintained at [github.com/ywqin829/gse276083-spatial-audit](https://github.com/ywqin829/gse276083-spatial-audit). A version-specific Zenodo DOI will identify the deposited analysis inputs and section-level results.

## Reproduce the audit

Extract analysis_code.zip into a folder named code, next to section_inputs.rds. From code, run the commands below. The primary run takes the November reference seed 42. The sensitivity run uses seeds 1 through 20. Tested with R 4.5.1, data.table 1.17.8, Matrix 1.7-3 and RANN 2.6.2.

```sh
Rscript --vanilla run_audit.R ../section_inputs.rds candidates.csv ../reproduced 42
Rscript --vanilla seed_sensitivity.R ../section_inputs.rds candidates.csv ../reproduced
```

The RDS holds ordered cell identifiers, section labels, broad and fine cell annotations, global coordinates in pixels, counts, segmented area, negative-probe counts per probe, neighbourhood covariates, cached reference distance and normalized expression for 19 unique genes. The candidate list specifies ten gene/cell-type pairs per disease and includes CXCL16 in both groups. The starting point for exact reproduction is this prepared input. The included export_inputs.R documents extraction from the November reannotated Seurat object, distance cache and neighbourhood table using Seurat 5.3.0. Reconstructing the original manual annotations from the GEO tables alone is outside this release.

## Archived outputs

The versioned archival record will contain the prepared analysis input and machine-readable output tables needed to reproduce and inspect the audit. Biological interpretation and result figures are intentionally not presented on this repository landing page.
