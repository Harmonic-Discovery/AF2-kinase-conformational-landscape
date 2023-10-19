# Investigating the conformational landscape of AlphaFold2-predicted protein kinase structures

All scripts used to generate the figures in the paper "Investigating the conformational landscape of AlphaFold2-predicted protein kinase structures" by Al-Masri et al. 2023 are included. 

![Alt text](outputs/figure.png?raw=true "Title")

All inputs needed are in the `inputs` directory and all outputs generated in `outputs`. The scripts include:

## AF2 Database Structure Conformations & Pocket Residue Fingerprints

#### 1. `AlphaFold_fingerprints.ipynb`
* Extracts the pocket residues of AF2 structures
* Generates KiSSim fingerprints for pocket residues
* Performs t-SNE on distance fingerprints

#### 2. `plot_conformations.ipynb`
* For both AF2 and PDB structures, the distribution of conformations by species & kinase group are plotted

#### 3. `plot_plddt.ipynb`
*  Plots average pLDDt for all pocket residues

#### Requirements 

  Python          # 3.6.8+
    biopython     # 1.77+
    pandas        # 1.4+
    numpy         # 1.22+
    plotly        # 5.8+
    sklearn       # 1.0.2+
    kissim        # 1.0.0

## Docking and MD simulations
#### 4. `plot_rmsd.py`, `plot_rmsf.py`, and `plot_MD_interaction_anaysis.py`
*  These script were used to create the plots to analyze the MD simulations

#### 5. `plot_docking_enrichment.py`
*  This file contains the functions used to compute and plot the enrichment from docking.

#### 6. `crossdocking_interaction_analysis.ipynb`
*  This notebook contains the functions used to create the plot the interactions analysis from docking. 

#### Requirements 

  Python          # 3.7.12+
    numpy         # 1.21.6
    pandas        # 1.4.2
    matplotlib    # 3.5.2
    seaborn       # 0.11.2
    MDAnalysis    # 2.1.0
    prolif        # 2.0.0
    biopython     # 1.79
    mdtraj        # 1.9.6
