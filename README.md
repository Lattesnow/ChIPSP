# ChIPSP: Spatial ChIP (ChIPSP) as a New Bioinformatics Tool to Characterize Spatial Gene Regulation

## Overview

`ChIPSP` is an R package implementing the **ChIPSP** workflow for integrating **ChIP-seq transcription factor binding** with **Hi-C chromatin interaction** data to identify spatially linked regulatory regions.

The current workflow contains five major steps:

1. Optional removal of chromosome X and chromosome Y interactions.
2. Merging Hi-C loop outputs across replicates and/or resolutions.
3. Spatial linking of ChIP-seq peaks to distal Hi-C loop anchors and ranking of ChIPSP interactions.
4. Gene annotation of conventional ChIP-seq peaks and ChIPSP regions.
5. Enrichr-based functional and pathway enrichment analysis.

The package is designed as a modular workflow so that individual components can also be run independently.

---

## Conceptual Framework

Hi-C identifies interacting chromatin regions, but the genomic bins used for Hi-C loop detection do not necessarily identify the exact nucleotide-level regulatory element responsible for the interaction.

ChIPSP integrates transcription factor ChIP-seq peaks with Hi-C loops to identify spatially linked regulatory regions.

For each Hi-C loop, ChIPSP determines whether a ChIP-seq peak overlaps either loop anchor. When a ChIP-seq peak overlaps one anchor, the interacting partner anchor is assigned as a spatially linked regulatory region.

The full interacting Hi-C anchor interval is therefore treated as a potential spatial regulatory region associated with the ChIP-seq transcription factor.

Spatial interactions are ranked using information from both:

- ChIP-seq peak strength (`pileup`)
- Hi-C interaction confidence (`FDR`)

The resulting ChIPSP regions can subsequently be used for gene annotation, pathway enrichment, transcription factor enrichment, genome-browser visualization, and other downstream analyses.

---

# Installation

## Reinstall the Current Development Version

If a previous version of `ChIPSP` is already installed, remove it first:

```r
remove.packages("ChIPSP")
```

Install the current development version directly from the `main` branch of GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "Lattesnow/ChIPSP@main",
  force = TRUE,
  upgrade = "never"
)
```

Load the package:

```r
library(ChIPSP)
```

The currently exported package functions can be checked with:

```r
getNamespaceExports("ChIPSP")
```

The current workflow includes the following major functions:

```r
removeXYChromosomes()
mergeHiCLoops()
readChIPFile()
chipSPLink()
ChIPSPannotation()
ChIPSPenrichment()
```

---

# Main Workflow

A typical complete ChIPSP analysis can be performed as follows:

```r
library(ChIPSP)

# ------------------------------------------------------------
# Step 1. Remove chromosome X and Y interactions
# ------------------------------------------------------------

hic_files_clean <- removeXYChromosomes(
  path = getwd()
)


# ------------------------------------------------------------
# Step 2. Merge Hi-C loop files
# ------------------------------------------------------------

hic_df <- mergeHiCLoops(
  hic_files_clean
)


# ------------------------------------------------------------
# Step 3. Perform ChIPSP spatial integration
# ------------------------------------------------------------

chipsp_results <- chipSPLink(
  hic_df = hic_df
)


# ------------------------------------------------------------
# Step 4. Annotate conventional ChIP and ChIPSP regions
# ------------------------------------------------------------

annotation_results <- ChIPSPannotation(
  chip_file = "r1881_hg19_test_peaks_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg19"
)


# ------------------------------------------------------------
# Step 5. Perform Enrichr pathway analysis
# ------------------------------------------------------------

enrich_results <- ChIPSPenrichment()
```

Each step is described in more detail below.

---

# Step 1: Remove Chromosome X and Chromosome Y

`removeXYChromosomes()` optionally removes genomic interactions located on chromosome X or chromosome Y.

This step can be useful when the analysis is intended to focus exclusively on autosomal chromatin interactions.

The function can automatically identify compatible Hi-C files in the specified directory.

Example:

```r
hic_files_clean <- removeXYChromosomes(
  path = getwd()
)
```

The function searches the working directory for compatible Hi-C files and generates corresponding files with chromosome X and chromosome Y interactions removed.

The returned object contains the paths to the filtered Hi-C files and can be passed directly to `mergeHiCLoops()`:

```r
hic_df <- mergeHiCLoops(
  hic_files_clean
)
```

This step is optional. If sex-chromosome interactions are biologically relevant to the study, users may skip this step and directly provide the original Hi-C files to `mergeHiCLoops()`.

---

# Step 2: Merge Hi-C Loop Files

`mergeHiCLoops()` combines Hi-C loop outputs from multiple files into a unified Hi-C loop table.

This allows ChIPSP to integrate loops obtained from:

- multiple biological replicates
- multiple Hi-C loop-calling resolutions
- multiple compatible Hi-C datasets

For example:

```r
hic_files_clean <- removeXYChromosomes(
  path = getwd()
)

hic_df <- mergeHiCLoops(
  hic_files_clean
)
```

The resulting `hic_df` is a merged Hi-C loop data frame that can be passed directly to `chipSPLink()`.

---

# Step 3: ChIPSP Spatial Integration

`chipSPLink()` performs the core ChIPSP analysis.

The function identifies ChIP-seq peaks overlapping Hi-C loop anchors and projects the ChIP-seq regulatory signal to the corresponding interacting anchor.

If a ChIP-seq peak overlaps anchor 1 of a Hi-C loop, anchor 2 is identified as the spatially linked region.

Likewise, if a ChIP-seq peak overlaps anchor 2, anchor 1 is identified as the spatially linked region.

Example:

```r
chipsp_results <- chipSPLink(
  hic_df = hic_df
)
```

When `chip_file` is not explicitly provided, `chipSPLink()` automatically searches for a compatible ChIP file in the working directory through `readChIPFile()`.

For example, a file such as:

```text
r1881_hg19_test_peaks_ChIP.bed
```

can be detected automatically.

Compatible ChIP-seq input formats include commonly used formats such as:

```text
BED
CSV
TSV
TXT
TAB
XLS
XLSX
```

The default output file is:

```text
ChIPSP_results.csv
```

The ChIPSP output contains genomic coordinates and ranking information including:

```text
chr
start
end
pileup
FDR
source_anchor
pileup_norm
fdr_norm
score
rank
```

The ChIPSP score integrates normalized ChIP-seq signal and normalized Hi-C interaction confidence.

Regions with stronger ChIP-seq signal and more confident Hi-C interactions receive higher ChIPSP ranking scores.

---

# Step 4: Gene Annotation

`ChIPSPannotation()` annotates both conventional ChIP-seq peaks and ChIPSP spatially linked regions to nearby genes.

The function uses genome-specific gene annotation resources and supports both human and mouse genomes.

Example for human hg19 data:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "r1881_hg19_test_peaks_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg19"
)
```

By default, two annotation files are generated:

```text
ChIP_anno_genes_upAdown_UCSC_Control.csv
ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv
```

The first file contains annotations for conventional ChIP-seq peaks.

The second file contains annotations for ChIPSP spatially linked regions.

---

## Species and Reference Genome Options

`ChIPSPannotation()` currently supports four species/genome-build combinations.

| Species | `species` | Supported `ref_genome` |
|---|---|---|
| Human | `"human"` | `"hg19"` |
| Human | `"human"` | `"hg38"` |
| Mouse | `"mouse"` | `"mm9"` |
| Mouse | `"mouse"` | `"mm10"` |

The species and reference genome must correspond to each other.

For example, these are valid:

```r
species = "human"
ref_genome = "hg19"
```

```r
species = "human"
ref_genome = "hg38"
```

```r
species = "mouse"
ref_genome = "mm9"
```

```r
species = "mouse"
ref_genome = "mm10"
```

Combinations such as:

```r
species = "human"
ref_genome = "mm10"
```

are not valid and will produce an error.

---

## Human hg19 Example

For ChIP-seq and Hi-C data aligned to the human hg19 genome:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg19"
)
```

This uses:

```text
TxDb.Hsapiens.UCSC.hg19.knownGene
org.Hs.eg.db
```

for gene annotation.

---

## Human hg38 Example

For human data aligned to hg38:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg38"
)
```

This uses:

```text
TxDb.Hsapiens.UCSC.hg38.knownGene
org.Hs.eg.db
```

---

## Mouse mm9 Example

For mouse data aligned to mm9:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "mouse",
  ref_genome = "mm9"
)
```

This uses:

```text
TxDb.Mmusculus.UCSC.mm9.knownGene
org.Mm.eg.db
```

---

## Mouse mm10 Example

For mouse data aligned to mm10:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "mouse",
  ref_genome = "mm10"
)
```

This uses:

```text
TxDb.Mmusculus.UCSC.mm10.knownGene
org.Mm.eg.db
```

---

## Annotation Distance

The default annotation window is:

```r
binding_bp = 5000
```

corresponding to a ±5 kb gene-association window.

For example:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg19",
  binding_bp = 5000
)
```

Users may modify this value according to the biological question and annotation strategy.

---

# Step 5: Enrichr Functional Enrichment Analysis

`ChIPSPenrichment()` performs Enrichr analysis using the gene annotations generated by `ChIPSPannotation()`.

Because the default filenames generated by `ChIPSPannotation()` match the default inputs expected by `ChIPSPenrichment()`, the enrichment analysis can normally be started simply with:

```r
enrich_results <- ChIPSPenrichment()
```

The function automatically reads:

```text
ChIP_anno_genes_upAdown_UCSC_Control.csv
ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv
```

and constructs five gene sets:

```text
ChIP
ChIPSP
ChIPSP_only
ChIP_only
Shared
```

where:

```text
ChIP
```

contains genes associated with conventional ChIP-seq peaks,

```text
ChIPSP
```

contains genes associated with ChIPSP spatial regions,

```text
ChIPSP_only
```

contains genes uniquely identified by ChIPSP,

```text
ChIP_only
```

contains genes uniquely identified by conventional ChIP-seq annotation,

and:

```text
Shared
```

contains genes identified by both methods.

---

## Default Enrichr Databases

The default enrichment workflow currently includes:

```r
c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023",
  "GO_Molecular_Function_2023",
  "KEGG_2021_Human",
  "Reactome_2022",
  "ChEA_2022",
  "MSigDB_Hallmark_2020"
)
```

The enrichment function generates individual Enrichr result tables as well as manuscript-oriented enrichment dotplots.

For the dotplots:

- pathways are ranked by nominal P value
- dot color represents `-log10(P)`
- dot size represents the number of overlapping genes
- the x-axis represents gene ratio
- gene ratio is calculated as overlapping genes divided by the total number of genes in the pathway

The default output directory is:

```text
Enrichr_dotplots
```

---

## Customizing Enrichment Databases

The Enrichr databases can be customized using the `databases` argument.

For example:

```r
enrich_results <- ChIPSPenrichment(
  databases = c(
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "ChEA_2022"
  )
)
```

The current default database set was designed primarily for the human ChIPSP workflow and includes:

```text
KEGG_2021_Human
```

For mouse analyses, users should select Enrichr databases appropriate for mouse genes and the intended biological analysis.

---

# Complete Human hg19 Example

A complete human hg19 workflow is:

```r
remove.packages("ChIPSP")

remotes::install_github(
  "Lattesnow/ChIPSP@main",
  force = TRUE,
  upgrade = "never"
)

library(ChIPSP)

getNamespaceExports("ChIPSP")


# Step 1
hic_files_clean <- removeXYChromosomes(
  path = getwd()
)


# Step 2
hic_df <- mergeHiCLoops(
  hic_files_clean
)


# Step 3
chipsp_results <- chipSPLink(
  hic_df = hic_df
)


# Step 4
annotation_results <- ChIPSPannotation(
  chip_file = "r1881_hg19_test_peaks_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg19"
)


# Step 5
enrich_results <- ChIPSPenrichment()
```

---

# Example for Human hg38

Only the annotation genome option needs to be changed when the input genomic coordinates are based on hg38:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "human",
  ref_genome = "hg38"
)
```

---

# Example for Mouse mm9

For mouse mm9 data:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "mouse",
  ref_genome = "mm9"
)
```

---

# Example for Mouse mm10

For mouse mm10 data:

```r
annotation_results <- ChIPSPannotation(
  chip_file = "example_ChIP.bed",
  chipsp_file = "ChIPSP_results.csv",
  species = "mouse",
  ref_genome = "mm10"
)
```

---

# Output Files

A standard ChIPSP analysis produces several classes of output.

The core ChIPSP spatial-integration output is:

```text
ChIPSP_results.csv
```

Gene annotation produces:

```text
ChIP_anno_genes_upAdown_UCSC_Control.csv
ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv
```

Enrichment analysis produces an:

```text
Enrichr_dotplots/
```

directory containing gene lists, pathway-enrichment tables, a combined enrichment result table, and enrichment dotplots.

---

# Input Requirements

All input files can be placed in the same working directory, although complete file paths may also be supplied.

Hi-C files should contain the genomic coordinates and FDR information required by the ChIPSP workflow.

ChIP-seq input files should contain genomic coordinates and ChIP-seq signal information, including the `pileup` signal used for ChIPSP ranking.

Users should ensure that the ChIP-seq and Hi-C datasets use the same reference genome.

For example, hg19 ChIP-seq coordinates should be combined with hg19 Hi-C coordinates and annotated using:

```r
species = "human"
ref_genome = "hg19"
```

Similarly, mm10 ChIP-seq and Hi-C data should be annotated using:

```r
species = "mouse"
ref_genome = "mm10"
```

Genome builds should not be mixed within the same ChIPSP analysis.

---

# Package Scope and Design

`ChIPSP` provides an integrated but modular framework for:

```text
Hi-C preprocessing
        ↓
Hi-C loop merging
        ↓
ChIP–Hi-C spatial integration
        ↓
ChIPSP ranking
        ↓
Gene annotation
        ↓
Functional enrichment
```

The core spatial-integration result remains independent of downstream interpretation. Users can therefore use the ChIPSP output with alternative genome annotations, enrichment platforms, transcription factor databases, genome browsers, or visualization tools if desired.

---

# Notes

ChIPSP performance depends directly on the quality of the underlying ChIP-seq and Hi-C datasets.

The ChIP-seq and Hi-C datasets should use compatible chromosome naming conventions and the same reference genome.

Hi-C loop detection is resolution-dependent, and the number and genomic size of detected loop anchors can influence the number of spatial regulatory regions identified by ChIPSP.

Users should therefore consider ChIP-seq quality, Hi-C sequencing depth, loop-calling resolution, and genome-build compatibility when interpreting ChIPSP results.

---

# Citation

If you use `ChIPSP` in your research, please cite the associated manuscript describing the ChIPSP method.

Citation details will be added once available.

---

# Contact

For questions, issues, or feature requests, please open an issue in the GitHub repository or contact:

`tianyi.zhou@childrens.harvard.edu`
