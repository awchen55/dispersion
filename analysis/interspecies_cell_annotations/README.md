# Cell Type Annotation

Cell type annotation pipeline for the diploid (human/chimp) and allotetraploid
single-cell RNA-seq datasets used throughout the dispersion divergence
project. Labels assigned here feed every downstream DD/DE and cis/trans
mechanism analysis in `../`.

## Files

| File | Description |
|------|-------------|
| `cell_annotations_diploid_lines.Rmd` | Annotates the integrated diploid human/chimp object with CelliD; also builds the shared reference UMAP and projects the allotetraploid onto it |
| `cell_annotations_allotetraploid.Rmd` | Builds the allotetraploid Seurat object from per-lane STARsolo outputs, annotates it with CelliD, and merges human/chimp allele counts into one combined object |

---

## Diploid pipeline (`cell_annotations_diploid_lines.Rmd`)

1. **Split into 5 chunks** — the integrated diploid object
   (`combined_1kFeatures_integrated.rds`, from Kenneth's SeuratData pipeline)
   is split into 5 roughly equal chunks so CelliD's MCA step doesn't run out
   of memory.
2. **Load reference signatures** — combines three marker-gene signature
   sources (FCA/Human Cell Atlas, Linnarsson brain atlas, Gilad lab cardiac
   signatures; see Annotation method below) plus a manually curated iPSC
   marker set.
3. **Run CelliD per chunk** — filters to protein-coding genes expressed in
   >5 cells, computes HVGs (excluding MT/ribosomal/cell-cycle genes), runs
   MCA (50 dims), and scores every cell against every reference signature
   via a hypergeometric enrichment test.
4. **Remap fine-grained labels** to the 20 broad cell types used in
   downstream analysis (see Cell type remapping below).
5. **Save full metadata** — joins cell type labels onto the full Seurat
   metadata table.
6. **QC plots** — cell counts by cell type/species, the >=1,000-cells subset
   used downstream, and per-donor composition within each cell type.
7. **Annotate Seurat objects** — writes the cell type labels back onto both
   the full combined object and each of the 5 chunks.
8. **Split by cell type** — one Seurat object per cell type (>=1,000 cells
   in both human and chimp).
9. **Export to AnnData** — counts, normalized data, metadata, gene/cell
   lists, and any reductions per cell type, as flat files for Python-side
   conversion.
10. **Verify exports** — checks exported `.mtx` cell counts against the
    annotation table.
11. **Reference UMAP panels** — fits one UMAP on the diploid data
    (human+chimp combined, on the existing cross-species-integrated PCA) and
    reuses it for the human-only and chimp-only panels; the allotetraploid
    is then *projected* into this same coordinate space via
    `FindTransferAnchors`/`MapQuery` (never independently refit), so all four
    panels are directly comparable.

## Allotetraploid pipeline (`cell_annotations_allotetraploid.Rmd`)

Each allotetraploid cell carries both a human and a chimp genome, and reads
are separated by allele of origin at alignment time (Kenneth's pipeline
produces three alignment flavors per lane: **humanized** — human allele
reads aligned to a humanized chimp genome, used for annotation — **human**,
and **chimp**).

1. **Load humanized data (8 lanes)** — reads STARsolo raw counts per lane,
   filters to cells with >1,000 features (`create_object()` helper), and
   writes counts to BPCells on-disk storage to keep memory manageable.
2. **QC** — nFeature/percent.mt/percent.rp/percent.malat1 violin plots by
   lane.
3. **Normalize** and find variable features for CelliD.
4. **Load reference signatures** — same three-source set as the diploid
   pipeline.
5. **Run CelliD** on the humanized data, split into 2 memory-sized chunks.
6. **Remap labels** to the same 20 broad cell types and save
   `tetraploid_cell_annotations.csv`.
7. **Load human and chimp allele objects** separately (8 lanes each, same
   BPCells-backed approach as step 1).
8. **Match barcodes** across all three objects (humanized, human, chimp) —
   only cells present in all three have usable allele-specific data.
9. **Transfer annotations** onto the human/chimp allele objects, add a
   species prefix to barcodes (`RenameCells`), and merge into one combined
   object.
10. **Filter genes and cells** — drop genes with zero expression in either
    species (can't support a cis/trans comparison) and cells whose barcode
    isn't shared between the human and chimp allele data; re-apply the
    >1,000-feature filter on the combined object.
11. **Normalize and save** the final filtered object
    (`tetraploid_combined_variablfeatures_filtered.rds`) plus its metadata.
12. **Compare cell counts** against the diploid data (diploid: >=1,000 cells
    in both species; tetraploid: >=300 cells — each tetraploid cell counted
    once via its human-allele row to avoid double counting).
13. **Split by cell type and export to AnnData**, same flat-file format as
    the diploid pipeline.

---

## Annotation method: CelliD

Both pipelines use [CelliD](https://bioconductor.org/packages/CelliD)
(Multiple Correspondence Analysis embedding + per-cell hypergeometric
enrichment test against reference gene signatures) rather than
cluster-then-label marker scoring.

**Reference signatures** (three sources combined):
1. **FCA** (Cao et al. / Human Cell Atlas, via Wenhe Lyu) — broad cell type
   markers, plus a manually added iPSC marker set.
2. **Linnarsson brain atlas** — replaces FCA's brain categories with more
   granular neuronal subtypes.
3. **CM / cardiac signatures** (Gilad lab, Reem) — cardiomyocyte/mesoderm
   lineage markers.

Cell cycle genes (`cc.genes.updated.2019`) and non-protein-coding genes are
excluded from every signature; MT/ribosomal genes are excluded from the HVG
set used for MCA.

**Assignment threshold**: a cell gets the label with the highest enrichment
score if `-log10(corrected p-value) >= 2` (i.e. corrected p <= 0.01);
otherwise it's marked `"unassigned"`.

## Cell type remapping

Both pipelines remap fine-grained reference categories to 20 broad cell
types before anything downstream sees them (e.g. Visceral neurons /
Chromaffin cells / Sympathoblasts -> ENS neurons; Bipolar cells /
Photoreceptor cells / Horizontal cells / Amacrine cells / Ganglion cells ->
Retinal neurons; Fibroblast / Mesangial cells -> Stromal cells; see the
`## Step 4/5: Remap...` chunks in each Rmd for the full mapping).

A cell type is only used in downstream DD/DE analysis if it has >=1,000
cells in **both** human and chimp diploid lines and >=300 cells in the
allotetraploid — this yields the 20 cell types referenced throughout the
project (`celltypes_keep`/`celltypes2` in the DD/DE analysis notebooks).

---

## Outputs

| File | Produced by | Description |
|------|-------------|-------------|
| `data/cell_annotations/diploid_cell_annotations.csv` | diploid | Barcode -> cell type |
| `data/cell_annotations/diploid_metadata.csv` | diploid | Full Seurat metadata + cell type |
| `data/cell_annotations/tetraploid_cell_annotations.csv` | allotetraploid | Barcode -> cell type |
| `data/cell_annotations/tetraploid_metadata.csv` | allotetraploid | Full metadata for the filtered combined object |
| `data/diploid_crossfilt_data/combined_1kFeatures_integrated_annotated.rds` | diploid | Annotated diploid Seurat object |
| `data/tetraploid/tetraploid_combined_variablfeatures_filtered.rds` | allotetraploid | Final filtered allotetraploid object |
| Per-cell-type `.rds` + AnnData flat-file exports | both | `data/diploid_crossfilt_data/celltype/`, `data/tetraploid/celltype/` and their `anndata_exports/` subfolders |
| Reference UMAP panels (`.pdf`) | diploid | `results/dispersion_paper/main/umaps/` |

---

## R environment

Both notebooks are run under the project's standard HPC environment:

```
module load R/4.4.1 gsl hdf5/1.12.0
export R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1
```

R version: **4.4.1** (2024-06-14, "Race for Your Life")

### Package versions

Confirmed via `packageVersion()` in the environment above (not from
documentation or memory):

| Package | Version | Source | Used for |
|---------|---------|--------|----------|
| Seurat | 5.3.0 | CRAN | Core single-cell object/analysis |
| SeuratObject | 5.2.0 | CRAN | Seurat's underlying object classes |
| CelliD | 1.14.0 | Bioconductor | MCA + hypergeometric cell type annotation |
| BPCells | 0.3.1 | GitHub (bnprks/BPCells) | On-disk sparse count storage for the allotetraploid pipeline |
| Matrix | 1.7.3 | CRAN | Sparse matrix I/O (`writeMM`, `readMM`, row/col sums) |
| dplyr | 1.1.4 | CRAN | Data wrangling |
| tidyr | 1.3.1 | CRAN | Reshaping (`pivot_wider`/`pivot_longer`, `complete`) |
| tibble | 3.3.0 | CRAN | `rownames_to_column`/`column_to_rownames` |
| stringr | 1.5.1 | CRAN | Barcode string matching (`str_starts`, `str_remove`, `str_trim`) |
| ggplot2 | 3.5.2 | CRAN | QC and UMAP plots |
| scales | 1.4.0 | CRAN | Plot axis formatting, color palette (`hue_pal`) |
| parallel | (base R 4.4.1) | base | `mclapply` for per-lane loading |

`CelliD` and `BPCells` are not on CRAN — install via
`BiocManager::install("CelliD")` and `BiocManager::install("BPCells")`
respectively.
