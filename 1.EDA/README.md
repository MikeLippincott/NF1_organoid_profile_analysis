# 1. Exploratory Data Analysis (EDA)

This module contains exploratory analyses of image-based morphological profiles from NF1 organoids imaged in 2D (max projection) and 3D.
All analyses are run across four profile types: 2D organoid, 2D single-cell, 3D organoid, and 3D single-cell, each with normal, feature selected (fs), aggregate (agg), and consensus data.

Some notebooks have not yet been utilized specifically for 2D vs. 3D organoid analysis.


## Shannon Entropy (Notebooks 12–13)

### Why compute it?
Not all morphological features carry meaningful variation as some are effectively constant across the dataset.
Even after feature selection
Computing per-feature entropy helps identify which features are informative and supports downstream feature selection.
It also lets us compare the overall information content of 2D vs. 3D profiles.

### What is Shannon entropy?
Shannon entropy quantifies the uncertainty or unpredictability in a distribution.
It captures both how spread out values are and how evenly they occur across different outcomes.
A distribution where all outcomes are equally likely has maximum entropy, while one concentrated on a single outcome has zero entropy.
Entropy is computed in bits (base-2 log).

For morphological features, this means:
a feature with high entropy varies widely across organoids and is likely more informative,
while a feature with low entropy is nearly constant across samples and unlikely to help distinguish conditions.

### How is it computed?
Each feature's values are binned into a histogram with 50 bins.
The count in each bin is divided by the total number of samples to create a probability distribution.
Shannon entropy is then computed from this distribution using `scipy.stats.entropy`, which sums `-p * log2(p)` for each bin with nonzero probability.
If values are spread evenly across all bins, entropy is high (maximum ~5.6 bits for 50 bins).
If values are concentrated in just a few bins, entropy is low (near 0 bits).

### Why 50 bins?
There is no single correct method for choosing the number of bins for entropy estimation — it depends on the data.
Common heuristics like the square root rule (`sqrt(n)`) provide a rough starting point but were designed for histogram visualization, not entropy computation.
For ~3000 single-cell profiles, the square root rule gives ~55 bins, which is close to our choice of 50.
50 bins is a practical, round number that provides enough resolution to capture distribution shape in both 2D and 3D profiles, which have different numbers of features.
With 3000 samples and 50 bins, each bin averages ~60 samples, which is sufficient for stable probability estimates.
The number of bins can be adjusted if needed, but the relative entropy rankings across features are expected to remain stable across reasonable bin counts.

### Graphs
Each plot is a violin + boxplot comparing the entropy distributions of 2D and 3D features.
The violin shows the full shape of the distribution, while the boxplot inside shows the median, quartiles, and whiskers.
Individual dots are outlier features with unusually high or low entropy.
Each dot represents a single feature, not a treatment or patient.