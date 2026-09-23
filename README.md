# prepTCGAdata

Prepares TCGA splice-junction data for the [SpliceRx](https://github.com/agiusp/SpliceRx)
app. Downloads a TCGA cohort's splice-junction counts and clinical metadata,
annotates junctions against a GENCODE reference, scores novel junctions, and
writes everything out as the `.rds` files SpliceRx's Data tab expects.

## Install

```r
# install.packages("devtools")  # if you don't have it
devtools::install_github("agiusp/prepTCGAdata")
```

## Use

```r
library(prepTCGAdata)
get_tcga_data("~/tcga-cohort-data/PAAD", "PAAD")  # any TCGA project code, e.g. "COAD", "KIRC", "STAD"
```

`get_tcga_data(ddir, cohort)` writes the following files to `ddir`, each
prefixed `TCGA_<cohort>_`:

| File | What it is |
|---|---|
| `junction_counts.rds` | Raw junction × sample count matrix, pulled from `recount3` |
| `junction_metadata.rds` | One row per junction — coordinates, strand, recount3/Snaptron annotation, and (after `annotate_sj()`) the overlapping GENCODE gene id/name |
| `novel_junction_counts_per_gene.rds` | Novel (unannotated) junction counts, rolled up per gene and sample |
| `novel_junction_RRS_scores.rds` | Each novel junction's relative read support against its gene's annotated junctions — `X / (X + median annotated-junction support)`, so bounded to `[0, 1)` |
| `novel_junction_counts_per_pathway.rds` | The per-gene novel-junction counts rolled up into xCell cell-type and (where available) ConsensusTME tumour-microenvironment gene-set scores |
| `sample_metadata.rds` | recount3's per-sample columns joined with GDC overall-survival fields and cBioPortal MSI status/scores |

The whole pipeline (junction annotation, novel-junction scoring, pathway
roll-up) runs against the GENCODE release 29 gene models bundled with the
package, since that's what `recount3` builds its TCGA junctions against — no
separate GENCODE download needed. `sample_metadata.rds` carries no
Good/Poor survivor labelling of its own; that stratification happens
interactively in SpliceRx's SJSurv tab, where the age bands and histology
grouping can be tuned to the cohort at hand.

See `?get_tcga_data` and the other exported functions (`annotate_sj`,
`count_novel_sjs`, `count_novel_sjs_per_pathway`, `novel_junction_RRS`) for
the individual steps.

## License

MIT
