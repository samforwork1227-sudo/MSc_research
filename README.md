# MSc_research — Membrane-aware structure preparation for single-pass receptors

A reproducible computational pipeline that converts continuous, single-chain
transmembrane-receptor models into **membrane-embedded, energy-minimised**
conformations with spatially separated extracellular (ECD), transmembrane (TM)
and intracellular (ICD) regions. The workflow couples **Boltz-2** structure
prediction, **UCSF ChimeraX** + **Biopython** topology editing, and **OpenMM**
membrane construction / minimisation, orchestrated with **Nextflow** on a
SLURM HPC cluster. FGFR2 is the primary test case; **FGFR1** is used to confirm
that the pipeline generalises to a closely related single-pass receptor.

> This repository accompanies the MSc dissertation *"Improving the structure
> prediction of proteins spanning lipid membranes for improved drug discovery"*
> (King's College London, MSc Applied Bioinformatics).

## Pipeline overview

```
FGFR2 / FGFR1 sequence
      |  Boltz-2 (co-folded with 20 POPC-like lipids)
      v
  full-length model (.cif)
      |  ChimeraX  (domain re-orientation, TM axis, torsions)
      v
  pre-MD model (.pdb)
      |  Biopython  (linker straightening: FGFR2 ECD-TM 359-377 -> 18.0 A; ICD 399-418 -> 25.0 A;
      |              FGFR1 ranges shifted to the 377-397 TM annotation)
      v
  clean_linkers.pdb
      |  OpenMM  (PDBFixer -> implicit min -> addMembrane POPC -> coarse -> fine min)
      v
  membrane-embedded, minimised model (.pdb)
```

## Repository layout

The repository is organised as a numbered lab notebook. Each folder is one
experiment; folder numbers match the **Experiment ID** column of
`experiment_master.xlsx`.

| Folder | Stage |
|--------|-------|
| `001`-`008_fgfr2_*` | Boltz-2 lipid-ligand exploration (ligand count, colouring, constraints) |
| `009_fgfr2_pipeline` | First Boltz-2 -> ChimeraX -> OpenMM pipeline |
| `010_fgfr2_20lipids` | 20-lipid model generated on the HPC |
| `011_fgfr2_straighten` | Biopython linker straightening |
| `012_fgfr2_finetune` | Linker fine-tuning / domain separation |
| `013_nextflow_pipeline` | Nextflow pipeline excluding membrane construction and system energy minimisation steps(`main.nf`, `nextflow.config`) |
| `014_membrane` | Explicit membrane construction in OpenMM (four sub-steps) |
| `015_fgfr1` | **Final consolidated Nextflow pipeline** (full pipeline applied to FGFR1) |
| `dissertation_figures.ipynb` | Plotting notebook for the pLDDT, PAE, energy-profile and runtime figures |
| `experiment_master.xlsx` | Master experiment log |

### Suggested target structure

```
MSc_research/
|-- README.md
|-- LICENSE
|-- .gitignore
|-- envs/                     # environment.yml / requirements.txt / install_*.sh
|-- pipeline/                 # final Nextflow pipeline (was 015_fgfr1)
|   |-- main.nf
|   `-- nextflow.config
|-- experiments/              # exploratory runs 001-015 (lab-notebook history)
|-- analysis/                 # dissertation_figures.ipynb and derived figures
|-- data/                     # small inputs (sequences, YAML) - NOT large outputs
`-- docs/                     # dissertation figures, energy table
```

## Requirements

Exact versions used in this study (see also the **Software used** sheet of
`experiment_master.xlsx` and Appendix Table 4 of the dissertation):

| Software | Version |
|----------|---------|
| Boltz-2 | 2.2.1 |
| UCSF ChimeraX | 1.11.1 locally; HPC container tag `maduprey/chimerax:1.5` |
| Biopython | 1.84 |
| OpenMM | 8.5.0 |
| Nextflow | 25.10.4 (build 11173) |
| Java / JDK | OpenJDK 17.0.16+8 |
| Python (Boltz-2 env) | 3.11.6 |
| CUDA | 12.2.1 |

Separate environments are used for each tool (create with the scripts in
`envs/`). Figure generation additionally requires **NumPy** and **Matplotlib**
(see `dissertation_figures.ipynb`).

> **Note on portability:** paths in `main.nf` (`params.project_root`,
> `params.*_venv`) are absolute HPC paths. Edit these, or override on the
> command line, before running elsewhere.

## Running the pipeline

```bash
cd 015_fgfr1
# edit paths in main.nf / nextflow.config for your system, then:
sbatch run_nextflow_fgfr1.sbatch          # HPC
# or, interactively:
nextflow run main.nf
```

For FGFR1, the same pipeline is used with the FGFR1 sequence, the TM-flanking
pocket contacts set to A:377 / A:397, and the ChimeraX stage and linker ranges
adjusted to the FGFR1 numbering (`015_fgfr1`).

The longer OpenMM membrane-construction / minimisation stages are submitted as
dependent SLURM jobs because of HPC wall-time limits.

## Figures

`dissertation_figures.ipynb` reads the pipeline outputs (Boltz-2 confidence
files, OpenMM stage energies, SLURM runtimes) and produces the per-residue
pLDDT profiles, PAE heatmaps, staged energy profiles and runtime charts used in
the dissertation, via NumPy and Matplotlib.

## Data and code availability

All scripts, configuration files and representative intermediate structures are
in this repository. Full Boltz-2 prediction outputs (`.npz`, `.pkl`) and SLURM
logs are large and are best excluded from version control (see `.gitignore`).

## Citation

If you use this workflow, please cite the accompanying MSc dissertation, Improve the structure prediction of proteins spanning lipid membranes for improved drug discovery, King’s College London (2026).
