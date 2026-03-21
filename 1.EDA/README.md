# 1. Exploratory Data Analysis (EDA)

This module contains exploratory analyses of image-based morphological profiles from NF1 organoids imaged in 2D (max projection) and 3D.
All analyses are run across four profile types: 2D organoid, 2D single-cell, 3D organoid, and 3D single-cell, each with normal, feature selected (fs), aggregate (agg), and consensus data.

Some notebooks have not yet been utilized specifically for 2D vs. 3D organoid analysis.


## Shannon Entropy (Notebooks 12–13)

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
Each feature's values are binned into a histogram with 50 bins, and the count in each bin is divided by the total number of samples to create a probability distribution. Shannon entropy is then computed from this distribution using `scipy.stats.entropy` it sums `-p * log2(p)` for each bin with nonzero probability. If values are spread evenly across all bins, entropy is high (maximum ~5.6 bits for 50 bins). If values are concentrated in just a few bins, entropy is low (near 0 bits).

### Why 50 bins?
50 bins gives enough resolution to capture the shape of each feature's distribution without making bins so small they become sparsely populated, which would artificially inflate entropy estimates. 50 bins is a reasonable middle ground where it can capture differences yet still minimize noise. However, the number of bins can be changed as see fit.

### Graphs
Each plot is a violin + boxplot comparing the entropy distributions of 2D and 3D features. The violin shows the full shape of the distribution, while the boxplot inside shows the median, quartiles, and whiskers. Individual dots are outlier features with unusually high or low entropy. Each dot represents a single feature, not a treatment or patient.