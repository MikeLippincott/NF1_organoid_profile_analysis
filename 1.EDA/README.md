# 1. Exploratory Data Analysis (EDA)

This module contains exploratory analyses of image-based morphological profiles from NF1 organoids imaged in 2D (max projection) and 3D.
All analyses are run across four profile types: 2D organoid-aggregated, 2D single-cell-aggregated, 3D organoid-aggregated, and 3D single-cell-aggregated.

All notebooks have a matching `.py`/`.R` script.
Some notebooks have not yet been utilized specifically for 2D vs. 3D organoid analysis.


## Shannon Entropy (Notebooks 12–13)
**Scripts:** `12.calculate_shannon_entropy.py`, `13.plot_shannon_entropy.R`

### What is Shannon entropy?
Shannon entropy measures how spread out a distribution is.
A feature with high entropy varies widely across organoids and is likely more informative.
A feature with low entropy is nearly constant across samples and unlikely to help distinguish conditions.
Entropy is computed in bits (base-2 log).

### Why compute it?
Not all morphological features carry meaningful variation — some are effectively constant across the dataset.
Computing per-feature entropy helps identify which features are informative and supports downstream feature selection.
It also lets us compare the overall information content of 2D vs. 3D profiles.

### How is it computed?
Each feature's values are binned into a histogram, and Shannon entropy is computed from the resulting probability distribution using `scipy.stats.entropy`.

### Why 50 bins?
50 bins gives enough resolution to capture the shape of each feature's distribution without making bins so small they become sparsely populated, which would artificially inflate entropy estimates. 50 bins is a reasonable middle ground where it can capture differences yet still minimize noise. However, the number of bins can be changed as see fit.

### Outputs
- `1.EDA/results/entropy/` — per-feature entropy values for each profile type (`.parquet`)
- `1.EDA/figures/entropy/` — violin + boxplots comparing entropy distributions between 2D and 3D


## Directory Structure

```
1.EDA/
├── notebooks/       # Jupyter notebooks
├── scripts/         # Python and R scripts
├── figures/         # Output figures organized by analysis type
└── results/         # Output data files organized by analysis type
```
