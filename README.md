# OutlierPCsSelection

This repository contains WDL workflows and R scripts for identifying expression outliers using principal component analysis (PCA). The approach computes z-scores across a grid of PC counts from RNA-seq gene expression data, enabling selection of the optimal number of PCs for outlier detection. It also includes a workflow for filtering rare variants from a VCF file within cis windows.

## Repository Structure

```
.
├── workflows/
│   ├── CallZscoresAcrossPCs.wdl   # Computes z-scores across a range of PC counts
│   └── FilterRareVariants.wdl     # Filters rare variants within cis windows from a VCF
├── scripts/
│   ├── ComputeFitAcrossPCs.R      # Main R script for PCA and z-score computation
│   └── utils/
│       ├── ComputeBetasPCs.R      # Utility: PCA, beta estimation, and z-score functions
│       └── ExpressionPreprocessing.R  # Utility: GCT loading and count normalization
└── envs/
    └── Dockerfile                 # Docker image definition
```

---

## Workflows

### `CallZscoresAcrossPCs.wdl`

Runs the `ComputeFitAcrossPCs.R` script to load a GCT-format gene expression count matrix, normalize expression values, compute principal components, and calculate z-scores across a grid of PC counts. The output is an RDS file containing z-scores for each evaluated number of PCs.

**Inputs:**

| Name      | Type | Description                             |
|-----------|------|-----------------------------------------|
| `GCTFile` | File | Gene expression count matrix in GCT format |
| `Memory`  | Int  | Memory to allocate (in GB)              |

**Outputs:**

| Name      | Type | Description                                        |
|-----------|------|----------------------------------------------------|
| `Zscores` | File | RDS file (`ZscoresAcrossPCs.rds`) containing a list of z-score matrices, one per evaluated PC count |

---

### `FilterRareVariants.wdl`

Filters a VCF file to retain only variants overlapping defined cis windows, using `bcftools view`.

**Inputs:**

| Name         | Type | Description                                            |
|--------------|------|--------------------------------------------------------|
| `VCF`        | File | Input VCF file                                         |
| `VCFIndex`   | File | Index file for the input VCF                           |
| `CisWindows` | File | BED-format file defining cis window regions            |
| `GnomadAFs`  | File | gnomAD allele frequency file (reserved for future filtering steps) |
| `Memory`     | Int  | Memory to allocate (in GB)                             |

**Outputs:**

| Name                       | Type | Description                                                                        |
|----------------------------|------|------------------------------------------------------------------------------------|
| `CisWindowRareVariants.gz` | File | Compressed VCF containing variants within the defined cis windows (rare variant allele frequency filtering using gnomAD is reserved for a future step) |

---

## Scripts

### `scripts/ComputeFitAcrossPCs.R`

Main entry-point script called by `CallZscoresAcrossPCs.wdl`. It orchestrates the following steps:

1. Loads and filters a GCT expression count matrix using `LoadGCTFile` and `FilterCountData`.
2. Normalizes counts to log-CPM using `NormalizeCountsCPMs` (via edgeR).
3. Computes principal components with `ComputePCs` (via PCAtools).
4. Estimates beta (regression weight) values per PC using `ComputePCBetaValues`.
5. Computes z-scores across a grid of PC counts (0 to 8900, step 500) using `ZScores_Kgrid`. In practice, the grid is automatically capped at the number of available PCs (i.e., values exceeding the sample count are excluded).
6. Saves the result as `ZscoresAcrossPCs.rds`.

**Command-line Arguments:**

| Argument    | Type   | Description                               |
|-------------|--------|-------------------------------------------|
| `--GCTFile` | string | Path to the input GCT gene expression file |

---

### `scripts/utils/ExpressionPreprocessing.R`

Provides utility functions for loading and preprocessing expression data:

- **`LoadGCTFile(PathGCT)`** – Reads a GCT file (skipping the two header lines) into a data table.
- **`FilterCountData(CountData, CountThresh=6, PropSamples=0.2)`** – Transposes the count matrix and retains only genes where at least 20% of samples have counts greater than 6.
- **`NormalizeCountsCPMs(CountData)`** – Calculates TMM normalization factors with edgeR and returns a log-CPM matrix (samples × genes).

---

### `scripts/utils/ComputeBetasPCs.R`

Provides utility functions for PCA-based outlier z-score computation:

- **`ComputePCs(SamplesByGenes)`** – Runs PCA (via PCAtools) on the expression matrix and returns the PCA result object.
- **`ComputePCBetaValues(PCScores, SamplesByGenes)`** – Estimates regression weights (beta values) of each PC on each gene using ordinary least squares.
- **`ZScores_Kgrid(Z, BetaPCs, SamplesByGenes, K_grid, gene_block=2000)`** – Computes residual z-scores for each number of PCs in `K_grid`. `Z` is the PC scores matrix (samples × PCs); `BetaPCs` contains the regression weights (PCs × genes). Genes are processed in blocks of `gene_block` for memory efficiency. Returns a named list of z-score matrices (samples × genes), one per value of K.

---

## Docker Environment

All workflows run inside the Docker image:

```
ghcr.io/aou-multiomics-analysis/OutlierPCsSelection:main
```

The image is built from `envs/Dockerfile` using `mambaorg/micromamba:1.5.3` as the base. The following packages are installed:

**R packages (via conda-forge/bioconda):**
- `tidyverse`, `data.table`, `optparse`, `arrow`, `RNOmni`, `R.utils`, `janitor`, `magrittr`
- `bioconductor-pcatools` (PCA)
- `bioconductor-edger` (normalization)

**Other tools:**
- `bcftools` (VCF filtering)

The utility R scripts (`ComputeBetasPCs.R`, `ExpressionPreprocessing.R`) are copied to `/opt/r/lib/` inside the container.

---

## Running the Workflows

The workflows are registered on [Dockstore](https://dockstore.org/) and can be launched from there or run directly using a WDL-compatible executor such as [Cromwell](https://cromwell.readthedocs.io/) or the [Terra](https://app.terra.bio/) platform.

### Example: Running `CallZscoresAcrossPCs` with Cromwell

1. Prepare a JSON inputs file (e.g., `inputs.json`):

```json
{
  "ComputeZscores.GCTFile": "/path/to/expression.gct",
  "ComputeZscores.Memory": 32
}
```

2. Run with Cromwell:

```bash
java -jar cromwell.jar run workflows/CallZscoresAcrossPCs.wdl --inputs inputs.json
```

### Example: Running `FilterRareVariants` with Cromwell

1. Prepare a JSON inputs file:

```json
{
  "FilterRareVariants.VCF": "/path/to/variants.vcf.gz",
  "FilterRareVariants.VCFIndex": "/path/to/variants.vcf.gz.tbi",
  "FilterRareVariants.CisWindows": "/path/to/cis_windows.bed",
  "FilterRareVariants.GnomadAFs": "/path/to/gnomad_afs.vcf.gz",
  "FilterRareVariants.Memory": 16
}
```

2. Run with Cromwell:

```bash
java -jar cromwell.jar run workflows/FilterRareVariants.wdl --inputs inputs.json
```
