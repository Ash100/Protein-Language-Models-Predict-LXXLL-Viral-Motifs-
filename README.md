# Protein Language Models Predict LXXLL Viral Motifs

**In silico prioritization of viral LXXLL-containing motifs as candidate mimics of host nuclear receptor coactivators, targeting the androgen (AR) and estrogen (ER) receptor ligand-binding domains**

[![License](https://img.shields.io/badge/license-TBD-lightgrey.svg)]()
[![Status](https://img.shields.io/badge/status-preprint%2Funder%20review-orange.svg)]()

---

## Overview

Nuclear receptor coactivators (e.g., NCOA1–4, CBP, MED1, PGC1α) engage the androgen receptor (AR) and
estrogen receptor (ER) via short, degenerate **LXXLL "NR box" motifs**, where L is leucine and X is any
amino acid. A number of viral proteins also contain LXXLL-like motifs, raising the hypothesis that certain
viral proteins could structurally mimic host coactivators and interact with steroid hormone receptor
ligand-binding domains (LBDs).

This repository contains the notebooks and supplementary data used to:

1. Retrieve candidate proteins for human-infecting viruses from NCBI and identify LXXLL-containing motifs.
2. Embed the resulting 25-residue motif windows, alongside 7 host coactivator control motifs
   (NCOA1–4, CBP, MED1, PGC1α), using the **ESM2** protein language model (`ESM2-t33_650M_UR50D`).
3. Project the embeddings with **UMAP** and benchmark multiple clustering algorithms to identify viral
   motifs occupying the same embedding-space neighborhood as known coactivators.
4. Rank candidates by cosine similarity to control motifs and prioritize a shortlist for structural modelling.
5. *(Downstream, in a separate notebook to be added)* Model candidate viral proteins against the AR/ER LBDs
   (Boltz-2 + AlphaFold2 relaxation) and rank predicted binding via MM-GBSA.

> **This is a hypothesis-generating, purely computational study.** Receptor binding, transcriptional
> modulation, agonistic/antagonistic activity, and biological relevance of any candidate in this repository
> are **unvalidated** and require experimental (in vitro/in vivo) confirmation. Reported binding scores are
> computational ranking estimates, not experimentally measured affinities.

---

## Repository contents

| File | Description |
|---|---|
| [`Viral_Protein_Motif_Search.ipynb`](./Viral_Protein_Motif_Search.ipynb) | Stage 1: retrieves viral protein sequences from NCBI and scans for LXXLL motifs |
| [`Embeddings_ESM2.ipynb`](./Embeddings_ESM2.ipynb) | Stage 2: ESM2 embedding of motif windows, UMAP projection, clustering benchmark, and control-vs-candidate similarity ranking |
| [`Sup_File_1.csv`](./Sup_File_1.csv) | Final curated motif dataset: 355 LXXLL-containing 25-residue peptide windows (`ID`, `Sequence`, `Virus Name`, `Protein Name`, `Motif Position`, `Motif Sequence`, `Extracted Motif`) |
| [`Supplementary_Table__1.csv`](./Supplementary_Table__1.csv) | Clustering algorithm benchmark (Silhouette, Davies-Bouldin, Calinski-Harabasz scores for KMeans, Agglomerative, Spectral, DBSCAN, HDBSCAN) |
| `README.md` | This file |
| `LICENSE` | *(to be added)* |

**Planned additions:** control motif sequence table (NCOA1–4, CBP, MED1, PGC1α with positions), structural
modelling / MM-GBSA notebook, final candidate ranking table, and representative structural models.

---

## Methodology summary

| Step | Method | Parameters used |
|---|---|---|
| Protein retrieval | NCBI Entrez, RefSeq | Max 20 records/virus (`retmax=20`); targeted subset, not full proteomes |
| Motif extraction | Regex `L[A-Z]{2}LL` | 25-residue window (10 aa upstream + LXXLL + 10 aa downstream); windows near a protein terminus may be shorter |
| Embedding | ESM2 (`ESM2-t33_650M_UR50D`) | Mean-pooled final transformer layer (layer 33) |
| Dimensionality reduction | UMAP | `n_neighbors=10`, `min_dist=0.1`, `metric='cosine'`, `random_state=42` |
| Clustering benchmark | KMeans, Agglomerative, **Spectral**, DBSCAN, HDBSCAN | `random_state=42` where applicable; Spectral Clustering selected as best-performing (highest Silhouette score) |
| Candidate ranking | Cosine similarity to control motifs, per-motif/per-cluster silhouette score | Top-10 nearest candidates exported per control motif |
| Structural modelling *(planned notebook)* | Boltz-2 + AlphaFold2 relaxation | Complexes vs. AR (residues 715–740) / ER (residues 358–376) LBDs; on-/off-target distance cutoff: *[TBD]* |
| Binding ranking *(planned notebook)* | MM-GBSA decomposition | Reported as a computational ranking score, not a measured affinity; exact settings: *[TBD]* |

Software package versions (`umap-learn`, `hdbscan`, `scikit-learn`) were not pinned at analysis time (installed
via unversioned `pip install` in Colab); please record the exact versions if reproducing this analysis, e.g.
via `pip show umap-learn hdbscan scikit-learn` in the environment used.

### Clustering benchmark results

| Method | Silhouette Score | Davies-Bouldin Index | Calinski-Harabasz Index |
|---|---|---|---|
| KMeans | 0.4705 | 0.6256 | 645.562 |
| Agglomerative | 0.4585 | 0.6399 | 624.3239 |
| **Spectral** | **0.4829** | 0.4755 | 494.3701 |
| DBSCAN | 0.2649 | 1.7288 | 397.5854 |
| HDBSCAN | 0.4086 | 0.3431 | 88.2246 |

Spectral Clustering was selected on the basis of the highest Silhouette score and is the clustering method
used for the candidate prioritization and structural modelling steps.

---

## Key findings (summary)

- A dataset of **355** LXXLL-containing 25-residue peptide windows across **212** unique viral proteins was
  curated and embedded alongside 7 host coactivator control motifs.
- Benchmarking five clustering algorithms identified **Spectral Clustering** as optimal (Silhouette = 0.4829).
- Candidate motifs were ranked by cosine similarity to control motifs in ESM2 embedding space to prioritize
  the shortlist carried forward to structural modelling.
- Structural modelling (planned notebook) identified 4 candidate on-target proteins against AR and 6 against
  ER, based on proximity to the established LXXLL-binding regions of each receptor's LBD.

> **Note:** the motif dataset (`Sup_File_1.csv`) spans **62** distinct virus names, though the manuscript
> text currently states 53 human-infecting viruses. This discrepancy is being reconciled with the
> corresponding authors before final publication — please treat the manuscript's virus count as provisional
> until updated.

---

## Requirements

- Python 3.9+, GPU runtime recommended for the ESM2 embedding step.
- `biopython`, `pandas`, `numpy`, `tqdm`, `requests`, `matplotlib`, `seaborn`, `scikit-learn`, `scipy`,
  `umap-learn`, `hdbscan`, `transformers`, `fair-esm`.

```bash
pip install biopython pandas numpy tqdm requests matplotlib seaborn scikit-learn scipy umap-learn hdbscan transformers fair-esm
```

## Usage

1. Run [`Viral_Protein_Motif_Search.ipynb`](./Viral_Protein_Motif_Search.ipynb) to retrieve viral protein
   sequences from NCBI and identify LXXLL motifs (or start directly from `Sup_File_1.csv`, the curated
   output of this step, for exact reproducibility independent of live NCBI re-querying).
2. Run [`Embeddings_ESM2.ipynb`](./Embeddings_ESM2.ipynb) to generate ESM2 embeddings, reduce dimensionality
   with UMAP, benchmark clustering methods, and rank candidates by similarity to host coactivator controls.
3. *(Planned)* Run the structural modelling notebook to generate AR/ER LBD complexes and MM-GBSA rankings.

> Two cells in `Embeddings_ESM2.ipynb` that generate synthetic/placeholder data for illustrative purposes are
> disabled by default (`RUN_DEMO_CELL = False`) and are not part of the reported analysis.

---

## Data availability

The motif dataset, control motif sequences (once added), clustering benchmark results, embedding/clustering
notebooks, and representative structural models (once added) are provided in this repository. A supplementary
table listing all included viruses, protein accessions, and protein names can be derived from
`Sup_File_1.csv`.

## Limitations

- The NCBI search retrieved a maximum of 20 protein records per virus (`retmax=20`) and therefore represents
  a curated subset rather than complete viral proteomes; some LXXLL-containing candidates may not have been
  captured.
- No sequence-scrambled or randomly sampled negative/background control motifs were included in this version
  of the analysis.
- Reported MM-GBSA values (once the structural modelling notebook is added) are computational ranking scores,
  not experimentally measured binding affinities.
- All findings are computational predictions and require experimental (in vitro/in vivo) validation.

## Citation

If you use this code or dataset, please cite:

> [Author list]. Protein Language Models Predict Candidate Viral Coactivator-Mimicking Motifs Targeting
> the Androgen and Estrogen Receptors: An *In Silico* Prioritization Study. *[Journal]*, [year]. DOI: [TBD]

## License

This project is licensed under [TBD — e.g., MIT License / CC-BY-4.0]. See [`LICENSE`](./LICENSE) for details.

## Contact

For questions about the data or code, please open an issue in this repository or contact the corresponding
authors listed in the manuscript.
