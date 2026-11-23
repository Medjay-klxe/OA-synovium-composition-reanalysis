# OA synovial transcriptomics: frozen computational release

Manuscript: **Cell composition adjustment and macrophage expression heterogeneity in osteoarthritis synovium: a cross-cohort re-analysis and semi-synthetic evaluation**.

This directory contains the analysis scripts, frozen processed inputs, gene definitions and derived tables used by the submitted manuscript. No repository URL or DOI has been invented. Create the authors' GitHub/Zenodo release from this directory after selecting a reuse license. A repository is not created by this packaging step.

## Datasets and comparisons

- GSE89408: primary bulk comparison, 25 non-OA and 19 OA samples after finite reconstruction fit correlation >= 0.6.
- GSE216651: 3 healthy and 3 OA donors; 7,492 macrophages; donor-level disease comparison and decomposition.
- GSE248453: 8 OA donors; 7 eligible for primary UMI-decile matching; all 8 in the separate quartile sensitivity analysis.
- GSE283080: 4 OA donors; within-donor accompanying expression.
- GSE152805: 3 OA synovial donors; within-donor accompanying expression and cluster + UMI, Scrublet sensitivity.
- GSE82107 and GSE143514: external bulk analyses.
- GSE246416, GSE233502, GSE294277 and E-MTAB-8322: audited additional datasets; their limitations/roles do not imply independent OA–healthy validation.

GEO accession pages follow https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152805 (replace the accession for other GEO studies). E-MTAB-8322 is an ArrayExpress study.

## Environment

The GSE152805 archived sensitivity result was rerun without changing its algorithm on 2026-11-23 under **Python 3.11.9 and Scrublet 0.2.3**. All 24 donor-by-gene contrasts were identical to the frozen table (maximum absolute difference 0). GPNMB mean = -0.2296 and MARCO mean = -0.4442 (rounded to four decimals); both were negative in 3/3 donors. Exact versions are in `environment/python_verified_20261123.json`. R and package versions, and archived R sessionInfo records, are in `environment/`.

Use a virtual environment and install `requirements.txt`. Arial is required to reproduce the final typography; fonts are not redistributed here. For R analyses install the packages listed in `environment/R_versions.txt` and make `Rscript` available, or set `RSCRIPT` to its executable. Set an UTF-8-capable locale when using Unicode filenames. The entry point sets `OA_REPO_ROOT` automatically. No computer-specific drive or home path is required by the release scripts.

## Quick start

From the repository root:

```sh
python reproduce.py figures
python reproduce.py matching
```

The figure stage uses the frozen tables and the original UMAP coordinates. It does not recluster cells or regenerate statistical findings. PNG previews produced by the relocated scripts were verified pixel-for-pixel against all six final manuscript figures and Fig. S1.

Other defined stages:

```sh
python reproduce.py external_bulk
python reproduce.py bulk
python reproduce.py candidates
python reproduce.py mixtures
python reproduce.py decomposition
python reproduce.py gse152805
python reproduce.py simulation
python reproduce.py gap
```

The simulation and gap stages run four batches of 25 replicates per scenario using the frozen seeds. They are computationally expensive and were **not** rerun merely for packaging. The old setup script supplies shared functions and simulated profiles; `paired_scores.R` is the submitted four-scenario evaluation. Do not substitute the legacy exploratory run's summary for the submitted paired-score results.

For `gse152805`, download only the synovial matrix files corresponding to GSM4626763, GSM4626764 and GSM4626765 from GSE152805 and put the `.matrix.mtx.gz` files in `data/external/GSE152805/滑膜原始矩阵/`. Their original filenames must begin with the respective GSM accession. Fixed QC cell indices, cluster assignments and target-expression matrices are supplied alongside that directory; the sensitivity uses the existing clusters. The full raw archive and large Seurat objects are not redistributed in this code package. Verification predictions and rerun logs are in `verification/GSE152805_rerun/`.

## Script/data mapping to the final manuscript

| Final item | Analysis code | Frozen source data |
|---|---|---|
| Fig. 1 | `scripts/figures/plot_fig12356.py` | Fixed study-design labels in the script |
| Fig. 2 | `scripts/analysis/bulk_gene_sets.R` | `data/components/bulk/` |
| Fig. 3 | `scripts/simulation/paired_scores.R`, `gap_sensitivity.R` | `data/simulation/`, `data/gap/` |
| Fig. 4; Table 2 | `scripts/analysis/three_component_decomposition.R`, `macrophage_scores_and_two_component.R`, `scripts/figures/plot_fig4.py` | `data/components/macrophages/`, `data/components/donor_audit/`, `data/umap/` |
| Fig. 5 | `scripts/analysis/candidate_models.R`, `controlled_mixtures.R` | `data/components/candidates/`, `data/components/mixtures/` |
| Fig. 6; Table 3 | `scripts/crosscohort/standardized_decile.py`, `quartile_GSE248453.py`, `scrublet_152805.py`, `filter_152805_doublets.py` | `data/crosscohort/`, `data/external/` |
| Additional file 1 | Gene sets, six-module membership, effective weights and matched-donor summaries | `definitions/Additional_file_1.xlsx` |
| Additional file 2: Fig. S1 | `scripts/figures/plot_s1.py` | Frozen UMAP, two-component decomposition and leave-one-donor-out tables |

`definitions/six_modules.csv` specifies 44 unique genes. In M3_antiox_buffer, NFE2L2 has direction +1 and the other 11 genes have direction -1. Each gene direction is applied once before equally weighting the six module means.

Some internal output table filenames retain the original analysis-stage prefixes (for example `Fig3A_preprocessing_coefficients.csv` under `candidates/` is now used in final Fig. 5A). The mapping above, not those internal historical prefixes, defines the current manuscript figure numbers. Obsolete figure renderers are not included.

## Statistical boundaries

Donors are biological replicates. OA-only cohorts support within-donor expression comparisons, not an independent OA–healthy disease comparison. The two- and three-component decompositions have different parameterizations. Simulated known composition denotes RNA contribution, not observed tissue cell counts. Decile matching is the primary Fig. 6 analysis; quartile and cluster + UMI/doublet-filtered results are separate sensitivities. Frozen numerical tables are supplied unchanged; recomputation should write to `results/` rather than overwrite `data/`.

## Verification and provenance

- `verification/GSE152805_rerun_comparison.csv`: exact rerun agreement.
- `verification/figure_pixel_comparison.csv`: relocated figure-script checks.
- `verification/script_provenance.csv`: source-script relative paths and original SHA-256 hashes.
- `SHA256SUMS.txt`: checksums of release files.

The figures and harmonized matching were rerun from this directory for release validation. All R scripts were syntax-checked. Full bulk, candidate, controlled-mixture, decomposition and simulation computations were not repeated at this packaging stage; their frozen processed inputs and derived outputs are supplied for audit.

External bulk analyses use `scripts/analysis/external_bulk.R` and `data/external_bulk/`. Frozen cohort-audit summaries are in `data/cohort_audit/`; `scripts/crosscohort/compare_cohorts.py` regenerates candidate-cell UMI and detection comparisons. The archived Chinese-language evidence audit is contextual documentation, not an additional pooled analysis.
