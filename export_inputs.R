# Export inputs used by the six-section analysis from the November annotation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
suppressPackageStartupMessages({library(Seurat); library(data.table)})
obj <- readRDS(args[1])
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
m <- obj[[]]
genes <- unique(fread(args[4])$gene)
sections <- c("18_72_TB", "Lung_TB_18-78", "Lung_TB_18-79",
              "18_76_NTM", "18_80_NTM", "18_81_NTM")
m <- m[m$tissue %in% sections, , drop = FALSE]
cells <- rownames(m)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
probes <- grep("^NegPrb", rownames(counts), value = TRUE)
stopifnot(length(probes) == 19L, nrow(m) == 337875L, !anyDuplicated(cells))
expr <- LayerData(obj, assay = "RNA", layer = "data")[genes, cells, drop = FALSE]
cache <- readRDS(args[2])
ids <- match(cells, rownames(cache))
stopifnot(!anyNA(ids))
x <- data.table(cell = cells, tissue = as.character(m$tissue),
                fine = as.character(m$celltype_fine), broad = as.character(m$celltype_broad),
                x = m$CenterX_global_px, y = m$CenterY_global_px,
                area = m$Area, ncount = m$nCount_RNA,
                neg_per_probe = as.numeric(Matrix::colSums(counts[probes, cells])) / length(probes),
                cached_distance = cache$d_Neural[ids])
z <- fread(args[3])
stopifnot(!anyDuplicated(z$cell))
covariates <- c("n_neighbours", "Epithelial", "Immune", "Stromal", "Unclassified")
x <- cbind(x, z[match(cells, z$cell), ..covariates])
saveRDS(list(cells = x, expression = expr, sections = sections), args[5])
