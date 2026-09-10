# Domainome reproduction: how to run

This folder reproduces the Domainome clinical-variant analysis (Figure 4a and 4b [all variants] of Beltran et al.). It uses variants pulled from the MultiomicsKG Neo4j database instead of the original UniProt/ClinVar/gnomAD downloads used by the original authors.

It does the same analysis three times, once per variant subset:

| Folder | Prefix | Variants used |
|---|---|---|
| `all_variants/` | `01` | All SNVs linked to Domainome proteins |
| `core_variants/` | `02` | Only SNVs that map to residues with `location = "core"` |
| `surface_variants/` | `03` | Only SNVs that map to residues with `location = "surface"` |

Each folder has the same 4-step pipeline:

```
0X-variants_extraction*.ipynb                    (Python, queries Neo4j)
        │  writes 0Xa-*.csv, 0Xb-*.tsv, 0Xc-*.tsv
        ▼
0X-05_clinical_variants_uniprot*.Rmd             (R) reads 0Xa → writes 0X-05a-*labels*.txt
0X-06_clinical_variants_clinvar*.Rmd             (R) reads 0Xb → writes 0X-06a-*labels*.txt
        ▼
0X-07_clinical_variants_classification*.Rmd      (R) reads 0X-05a, 0X-06a, 0Xc → figures in output_files/
```

Run the steps in this order. Only after running step 07 the figure files will be available in output_files/

---

## 1. One-time setup

### 1.1 Download the original Domainome analysis files

The notebooks and Rmd files read the published Domainome `analysis_files` folder. Download it from Zenodo ([record 14356805](https://zenodo.org/records/14356805)):

```bash
mkdir -p ~/Downloads/domainome_zenodo && cd ~/Downloads/domainome_zenodo
curl -L -o files_to_reproduce_analyses.tar.gz \
  "https://zenodo.org/records/14356805/files/files_to_reproduce_analyses.tar.gz?download=1"
tar -xzf files_to_reproduce_analyses.tar.gz
```

The archive keeps the authors' original directory structure, so `analysis_files` is nested several levels deep (under `…/8_zenodo_submission/files_to_reproduce_analyses/`). Find it with:

```bash
find ~/Downloads/domainome_zenodo -type d -name analysis_files
```

Copy that folder into your Domainome folder, which you will point `DOMAINOME_FOLDER` at in §1.2. The result must be `$DOMAINOME_FOLDER/analysis_files/`:

```bash
mkdir -p "$DATASTORE/domainome"
cp -R "<path printed by find>" "$DATASTORE/domainome/"
```

The pipeline uses these files from `$DOMAINOME_FOLDER/analysis_files/`:

| Used by | Files |
|---|---|
| Notebooks | `clinical_variants_uniprot_proteinAPI_domainomegenes.txt`, `gnomad_v4_variants_domainome.txt`, `variant_summary_jan24_SNVs_coding.txt` |
| Rmd files | `domain_QC_summary_reproducibility_ranked.txt`, `mutated_domainome_merged_filtered.txt`, `mutated_domainome_merged_filtered_all_VEPs_final.txt`, `geneid_to_uniprotID_human.txt`, `domains_with_clinvar_mismatches.csv`, `gnomad_uniprotID_ensID_geneID`, `pfam_ID_description_table.txt`, `modes_of_inheritance_mech.csv`, `PFAM_ID_to_SCOP_class.tsv` |

### 1.2 Configure paths in `.env`

The notebooks and Rmd files read them from two environment variables, which are defined in the repository's `.env` file:

```bash
DATASTORE="/path/to/your/datastore"
DOMAINOME_FOLDER="$DATASTORE/domainome"     # folder that contains analysis_files/
MULTIOMICSKG_FOLDER=$(pwd)                  # repository root
```

| Variable | Used by | Purpose |
|---|---|---|
| `DOMAINOME_FOLDER` | notebooks and Rmd files | Location of the Domainome `analysis_files/` |
| `MULTIOMICSKG_FOLDER` | Rmd files | Repository root. The Rmd files derive `experiment_dir` from it, as `$MULTIOMICSKG_FOLDER/experiments/domainome_reproduction/<subset>` |
| `NEO4J_PWD`, `NEO4J_DATABASE` | notebooks | Neo4j credentials. Optional; the defaults are `multiomicskg` / `neo4j` |

**How to load the variables.** `.env` uses shell syntax (`$DATASTORE`, `$(pwd)`). VS Code's automatic `.env` loading and R's `.Renviron` cannot expand that syntax. Instead, source the file **from the repository root** in a terminal:

```bash
cd /path/to/MultiomicsKG
set -a; source .env; set +a                 # export every variable in .env
echo "$DOMAINOME_FOLDER" "$MULTIOMICSKG_FOLDER"   # check both are set
ls "$DOMAINOME_FOLDER/analysis_files" | head
```

Then start every tool (VS Code, RStudio, `jupyter`, `Rscript`) **from that same terminal**, so the tool inherits the variables. §2 and §3 show how.

### 1.3 Neo4j with MultiomicsKG loaded
- Neo4j must be running at `bolt://localhost:7687` with the MultiomicsKG loaded.

### 1.4 Python environment
Install the project's Poetry environment from the repository root. It provides `pandas`, `py2neo` and `ipykernel`:

```bash
poetry install
poetry env info --path      # note this path; it is the Jupyter kernel to select
```

### 1.5 R environment
You need R (≥ 4.x), pandoc, and these packages:

```r
install.packages(c("rmarkdown", "knitr", "ggplot2", "data.table", "stringr", "viridis",
                   "ggpubr", "ggrepel", "pROC", "PRROC", "mltools", "reshape2"))
```

RStudio ships with its own pandoc. To render from a terminal you need pandoc on your `PATH` (for example `brew install pandoc`).

---

## 2. Step 1: run the Jupyter notebook (variant extraction)

Do this in the folder of the subset you want (`all_variants`, `core_variants` or `surface_variants`).

**Important:** each notebook writes its outputs to its **current working directory**, using relative file names. The notebook must therefore run from inside its own folder, because the Rmd files look for the outputs there.

### Option A: VS Code (recommended)
1. Quit VS Code if it is running. An already-running instance does not pick up new environment variables.
2. From the repository root, load `.env` and open VS Code:
   ```bash
   set -a; source .env; set +a
   code .
   ```
3. Make sure Neo4j is running (see §1.3).
4. Open, for example, `experiments/domainome_reproduction/all_variants/01-variants_extraction.ipynb`.
5. Click **Select Kernel** (top right) → **Python Environments** → choose the Poetry env from §1.4.
6. Click **Run All**.

If you can't restart VS Code, add a cell at the top of the notebook and run it first:

```python
%env DOMAINOME_FOLDER=/path/to/your/datastore/domainome
```

### Option B: command line
`nbconvert` is not a project dependency, so you may need to add it to the environment first:

```bash
cd /path/to/MultiomicsKG
set -a; source .env; set +a
poetry run pip install nbconvert          # once
cd experiments/domainome_reproduction/all_variants
poetry run jupyter nbconvert --to notebook --execute --inplace 01-variants_extraction.ipynb
```

### Expected outputs (in the same folder)

| Subset | Files produced |
|---|---|
| all | `01a-clinical_variants_uniprot_proteinAPI_domainomegenes_from_multiomicskb.csv`, `01b-variant_summary_from_multiomicskb.tsv`, `01c-gnomad_v2_variants_from_multiomicskb.tsv` |
| core | `02a-…_from_multiomicskb_core.csv`, `02b-variant_summary_from_multiomicskb_core.tsv`, `02c-gnomad_v2_variants_from_multiomicskb_core.tsv` |
| surface | `03a-…_from_multiomicskb_surface.csv`, `03b-variant_summary_from_multiomicskb_surface.tsv`, `03c-gnomad_v2_variants_from_multiomicskb_surface.tsv` |

- `a` = UniProt-style clinical variants
- `b` = ClinVar-style variant summary
- `c` = gnomAD variants

All three files are tab-separated with no header, even the `.csv` one.

---

## 3. Steps 2–4: knit the Rmd files

### 3.1 Before you knit
1. Make sure `DOMAINOME_FOLDER` and `MULTIOMICSKG_FOLDER` are set in the R session (see below).
2. Make sure the `output_files/` sub-folder exists next to the Rmd files. It is already in the repository. If it is missing:
   ```bash
   mkdir -p experiments/domainome_reproduction/all_variants/output_files
   ```
3. Check that the Step 1 outputs (`0Xa`, `0Xb`, `0Xc`) are present in the folder.

### 3.2 Order
| # | all_variants | core_variants | surface_variants | Produces |
|---|---|---|---|---|
| 2 | `01-05_clinical_variants_uniprot.Rmd` | `02-05_…_uniprot_core.Rmd` | `03-05_…_uniprot_surface.Rmd` | `0X-05a-clinical_variants_measured_uniprot_proteinAPI_labels*.txt` |
| 3 | `01-06_clinical_variants_clinvar.Rmd` | `02-06_…_clinvar_core.Rmd` | `03-06_…_clinvar_surface.Rmd` | `0X-06a-clinical_variants_measured_clinvar_labels*.txt` |
| 4 | `01-07_clinical_variants_classification.Rmd` | `02-07_…_classification_core.Rmd` | `03-07_…_classification_surface.Rmd` | Figure PDFs in `output_files/` |

Each knit also writes a `.html` report next to its Rmd file.

### Option A: RStudio
1. Quit RStudio if it is running.
2. Start RStudio from a terminal with `.env` loaded. Run the RStudio binary directly; on macOS, an app started from the Dock or with `open -a` does not inherit the shell's variables:
   ```bash
   cd /path/to/MultiomicsKG
   set -a; source .env; set +a
   /Applications/RStudio.app/Contents/MacOS/RStudio experiments/domainome_reproduction/domainome_reproduction.Rproj &
   ```
   Alternatively, start RStudio normally and set the variables in the console before you knit. Knit inherits them from the session:
   ```r
   Sys.setenv(DOMAINOME_FOLDER    = "/path/to/your/datastore/domainome",
              MULTIOMICSKG_FOLDER = "/path/to/MultiomicsKG")
   ```
3. Check that `Sys.getenv("DOMAINOME_FOLDER")` in the console prints the expected path.
4. In the Files pane, go into `all_variants/` (or `core_variants/` / `surface_variants/`).
5. Open the `0X-05…Rmd` file and click **Knit** (or press ⇧⌘K). Wait for the HTML preview.
6. Do the same for `0X-06…Rmd`, then `0X-07…Rmd`.

> Use **Knit**, not "Run All" in the console. Each Rmd calls `setwd(base_dir)`, which is `$DOMAINOME_FOLDER`, in its first chunk. When you knit, knitr resets the working directory to the Rmd's folder after that chunk, so `ggsave("output_files/…")` writes to this folder's `output_files/`. If you run the chunks interactively, the working directory stays at `$DOMAINOME_FOLDER`, and the figures are written to `$DOMAINOME_FOLDER/output_files/`.

### Option B: terminal
Render each file from inside its folder:

```bash
cd /path/to/MultiomicsKG
set -a; source .env; set +a
cd experiments/domainome_reproduction/all_variants
Rscript -e 'rmarkdown::render("01-05_clinical_variants_uniprot.Rmd")'
Rscript -e 'rmarkdown::render("01-06_clinical_variants_clinvar.Rmd")'
Rscript -e 'rmarkdown::render("01-07_clinical_variants_classification.Rmd")'
```

For the other subsets, change the folder and the prefixes, for example `02-05_clinical_variants_uniprot_core.Rmd`.

---

## 4. Results

The main reproduced figures are in each subset's `output_files/` folder:

- `Figure4a_aPCA_vs_pathogenicity_density.pdf`
- `Figure_4b_destabilization_vs_pathogenicity_proportions_allvariants.pdf`

The file names are the same in every subset folder, including `core_variants` and `surface_variants`.

The `0X-07` step also writes the other Figure 4 / ED Figure 4 panels (ROC curves, per-family plots, MECP2/CRX examples, ESM1v variance-explained plots, and so on) to the same folder.

> ⚠️ **Side effect:** `0X-07` also writes `number_of_pathogenic_by_family.txt` and `clinical_variants_with_fitnessdata.txt` into `$DOMAINOME_FOLDER/analysis_files/`. This **overwrites the downloaded Domainome files of the same name**, and every subset overwrites the previous subset's copy. If you need the originals, keep the extracted Zenodo archive, or point those two `write.table` calls (around lines 52–56) at `experiment_dir`.

---

## 5. Quick checklist (one subset)

```text
[ ] Zenodo analysis_files copied to $DOMAINOME_FOLDER/analysis_files
[ ] .env sourced from the repository root (DOMAINOME_FOLDER, MULTIOMICSKG_FOLDER set)
[ ] Neo4j running on bolt://localhost:7687, NEO4J_PWD set if not default
[ ] Run 0X-variants_extraction*.ipynb from inside its folder → 0Xa, 0Xb, 0Xc files created
[ ] Knit 0X-05 → 0X-05a labels file
[ ] Knit 0X-06 → 0X-06a labels file
[ ] Knit 0X-07 → PDFs in output_files/
```

## 6. Troubleshooting

| Symptom | Likely cause |
|---|---|
| `FileNotFoundError: 'None/analysis_files/…'` in the notebook | `DOMAINOME_FOLDER` was not set when VS Code/Jupyter started. See §2 |
| `cannot change working directory` at `setwd(base_dir)` in an Rmd | `DOMAINOME_FOLDER` is empty in the R session. See §3 Option A |
| `File '/experiments/domainome_reproduction/…' does not exist` | `MULTIOMICSKG_FOLDER` is empty in the R session |
| `File '…/0Xa-…csv' does not exist` in an Rmd | The notebook was not run, or was run from another working directory |
| `ServiceUnavailable` / connection refused in notebook cell 2 | Neo4j is not running, or the bolt port is not 7687 |
| `AuthError` | `NEO4J_PWD` is wrong, or was not set in the environment that launched Jupyter/VS Code |
| Empty DataFrames for core/surface | The DB has no `maps_res`/`residue.location` data loaded |
| `cannot open file 'output_files/…pdf'` | `output_files/` is missing, or the chunks were run interactively instead of knitted |
| `pandoc version 1.12.3 or higher is required` | Rendering from a terminal without pandoc on your `PATH` |
