"""Background text for every module tab, section, and sub-panel of the viewer.

Each entry says what the data are and how we derive them. Text is plain
markdown shown in a collapsed expander.
"""

import streamlit as st

PROFILE_GLOSSARY = """
**Profile types used throughout**
- **sc** (single-cell): one row per cell or nucleus the image analysis finds.
- **organoid**: one row per organoid the image analysis finds.
- **fs** (feature-selected): standardized measurements that remain after
  `pycytominer` drops the least useful ones.
- **agg** / **consensus**: summaries of the feature-selected profiles, one row per
  replicate or per treatment (not per cell).
- **2D** profiles come from a single image slice (either the sharpest point in
  the stack or the middle slice); **3D** profiles come from the full stack of
  depth images.
- We exclude `NF0037_T1_CQ1` from every EDA result (it is a separate analysis).
"""

MODULE_BACKGROUND = {
    "Overview": (
        "Experiment design (not computed): which drug went into which well, and "
        "which patient tumor samples ran on each platemap."
    ),
    "EDA": (
        "Exploratory analysis of the image-based shape and appearance measurements of NF1 "
        "patient tumor organoids treated with drugs at one or more doses (DMSO is the "
        "control). Every table here comes precomputed; the viewer never recomputes "
        "anything.\n" + PROFILE_GLOSSARY
    ),
    "Viability prediction": (
        "Can organoid shape and appearance predict cell viability? An Elastic Net "
        "regression (a statistical model that predicts an outcome from many measurements "
        "at once, while automatically ignoring the least useful ones) learns from "
        "image-based profiles to predict per-well viability. The "
        "target is the **per-patient min-max scaled viability**: we scale it so that no "
        "information leaks from the test data (the min and max come only from the "
        "training rows of each split). The model puts features on the same scale, "
        "and `ElasticNetCV` automatically chooses its two tuning settings "
        "(`alpha` and `l1_ratio`) by trying many combinations.\n\n"
        "We split the data into training and test sets three ways:\n"
        "- **LOPO**: leave one patient out (train on the other patients).\n"
        "- **LOTO**: leave one treatment out.\n"
        "- **Random split**: repeated random train/test splits.\n\n"
        "We also fit each model with **shuffled** (randomly scrambled) labels "
        "(`shuffle_status`) as a baseline for comparison. A useful model should do "
        "better than this scrambled-label version."
    ),
    "Linear modeling": (
        "Per-feature linear models that ask which shape and appearance measurements change "
        "with treatment. We fit one linear regression model per "
        "**(patient, treatment/dose, feature)** on 3D standardized profiles:\n\n"
        "`feature ~ treatment + cell_count + organoid_count + cell_per_organoid_count`\n\n"
        "The **technical** model set adds location-related factors: well distance "
        "from the plate center, cell x/y/z position and depth in the image stack. "
        "The **original** model set leaves those out. We fit models on "
        "organoid and single-cell profiles, and on well-level averages "
        "(median, `*_agg`). We adjust the p-values to account for testing many "
        "features at once (`pvalue_fdr`). "
        "The `coefficient` is the size of the effect, measured in standard "
        "deviations of the feature; for `treatment` it is the difference from DMSO."
    ),
    "Treatment vs DMSO": (
        "We compare treatment vs DMSO control per patient. We measure morphology "
        "differences as the **mean absolute difference from DMSO** across features "
        "(using the normalized scale). We compute viability as "
        "`log2(treatment % / DMSO %)`per patient. MEK signatures pool the "
        "four MEK inhibitors (Binimetinib, Mirdametinib, Selumetinib, Trametinib) "
        "at 1 and 10 µM.\n" + PROFILE_GLOSSARY
    ),
    "QC checks": (
        "Do the Treatment vs DMSO comparisons hold up against the plate itself, "
        "rather than the biology? We check that where a well sits on the plate "
        "doesn't drive the morphology comparisons, and that each well's profile "
        "still looks like DMSO when it should."
    ),
}

SECTION_BACKGROUND = {
    # ---- 0.Overview ----
    "Platemap": (
        "The well -> treatment/dose layout of each named platemap. The screen "
        "used two layouts; pick one to see its grid and which patient samples "
        "ran on it."
    ),
    "Drugs": (
        "Every drug (excluding the DMSO control), with its dose(s) and "
        "mechanism of action (how it works in the cell). The chart shows each "
        "drug's FDA approval status, colored by mechanism of action."
    ),
    "Patients & tumor manifestations": (
        "Every patient tumor sample screened, which platemap it ran on, and "
        "its tumor manifestation (cNF, pNF, MPNST or Other). The chart shows "
        "how many patients fall into each manifestation."
    ),
    # ---- 1.EDA ----
    "UMAP": (
        "2D UMAP layouts of single-cell shape and appearance measurements. UMAP is a "
        "technique that takes many measurements per cell and arranges them on a 2D "
        "plot so that similar cells end up near each other. "
        "**Pooled** layouts use all patients together, "
        "per projection (2D max, 2D middle, 3D) and profile variant. We build "
        "**Per-patient** layouts separately for each patient, so UMAP1/UMAP2 "
        "are comparable only *within* a patient, never across patients. "
    ),
    "PCA": (
        "PCA layouts of the pooled (all-patient) feature-selected, aggregated, "
        "and consensus profiles for 2D and 3D. PCA "
        "(Principal Component Analysis) finds the few main patterns that explain "
        "most of the variation across many measurements. Columns "
        "`PC0…PCn` are a sample's score on each of these patterns; the caption and "
        "scree plot show how much of the total variation each pattern captures. "
        "Use the X/Y component selectors to move between patterns."
    ),
    "Correlation heatmaps": (
        "Sample-by-sample **Pearson correlation** matrices of profiles "
        "(a standard way of measuring how similar two sets of numbers are). A high "
        "value means two samples have similar shape and appearance measurements."
    ),
    "Cell counts": (
        "Number of cells segmented per organoid, patient, and treatment, counted "
        "from the single-cell profiles for each profile type and split by 2D and 3D. "
        "Use it to check whether a treatment changes "
        "cell number; the linear models also account for this factor."
    ),
    "Neighbors": (
        "How closely nuclei and organoids pack together. "
        "**2D**: `Nuclei_Neighbors_*` (nearest and second-nearest distance, number "
        "of neighbors, percent touching, angle between neighbors) and "
        "`Organoid_Neighbors_NumberOfNeighbors_Adjacent`. **3D** has no equivalent "
        "columns, so it reports different measures: neighbors adjacent in the local "
        "shell, distance from the organoid center and exterior, and nucleus volume. "
        "We approximate 3D organoid neighbor counts: two organoids in the same "
        "well/field of view count as adjacent if their enclosing spheres "
        "(center point plus average outer radius, in µm) touch. Do not compare 2D "
        "and 3D columns directly."
    ),
    "Intensity": (
        "Mean and median fluorescence intensity per channel and compartment. "
        "The whole-organoid compartment comes "
        "from organoid tables; cell and nucleus compartments come from "
        "single-cell tables. We standardize `value` **per patient** (within patient × "
        "compartment × channel × statistic, across all treatments), so panels "
        "share one scale and you read each treatment relative to that patient."
    ),
    "Count vs viability": (
        "Mean cells per organoid (3D) joined with measured viability from the "
        "platemap, per patient × treatment × dose. We drop and log patients "
        "present in only one of the profiles or the platemap."
    ),
    "Consensus heatmaps": (
        "The same sample-by-sample correlation matrices as **Correlation heatmaps**, "
        "preset to the **consensus** profile type (one profile per replicate or "
        "treatment, summarized from single cells). "
        "Pick another profile type or normalization in the controls to compare."
    ),
    "Correlation vs viability": (
        "Each point is a **pair** of samples. **x**: how similar their shape and "
        "appearance measurements are (correlation). **y**: how much their "
        "viability differs (a magnitude, not signed -- high means the pair's "
        "viabilities are very different, not that either one is high or low). "
        "Most pairs are an unremarkable **Typical pair**. The labelled corners are "
        "rare pairs past both cutoffs: blue means morphology and viability agree "
        "(similar/similar or different/different, as expected if morphology "
        "tracks viability); red means they disagree (similar morphology but "
        "different viability, or the reverse) -- the surprising cases."
    ),
    # ---- 3.viability_prediction_models ----
    "Model performance": (
        "Elastic Net performance for each profile type, split method, and shuffle "
        "status. **Fold metrics** are one row per train/test fold; **summary "
        "metrics** average over folds. Metrics are R² (the share of variation the "
        "model explains; can be negative if it does worse than just "
        "guessing the average), and MSE and MAE (two measures of average prediction "
        "error, where lower is better). Compare each real model to its "
        "shuffled-label control."
    ),
    "Predicted vs actual": (
        "Held-out predictions from the models: each row is a sample with its "
        "`Actual_Viability` (per-patient min-max scaled) and `Predicted_Viability` "
        "from a model that did not see it in training. Points on the diagonal "
        "are perfect predictions. The viewer doesn't load the wide feature columns here."
    ),
    "Feature importances": (
        "Elastic Net coefficients per feature, split by profile type, split method, "
        "shuffle status, and image mode. The model puts features on the same scale before "
        "fitting, so coefficient size is comparable across features. Zero means "
        "the model dropped the feature as unhelpful. The "
        "sign shows the direction (a positive coefficient means higher viability)."
    ),
    # ---- 4.linear_modeling ----
    "Volcano (effects vs significance)": (
        "One point per model term: `coefficient` (effect size) against a "
        "transformed version of the adjusted p-value (`-log10(pvalue_fdr)`) "
        "where taller points are more statistically significant. A point is "
        "*significant* when `pvalue_fdr < 0.05`. "
        "Split by `term` to separate treatment from the other factors. Pick a "
        "model result file (organoid or sc; original or technical; agg or not) first."
    ),
    "Effect sizes": (
        "Distribution of model `coefficient` values by treatment or therapeutic "
        "category. For the `treatment` term this is the change versus DMSO, "
        "measured in standard deviations of that feature. The other terms "
        "(counts, position) are per-unit slopes and are not on the same scale."
    ),
    "Model fit": (
        "How well each model fits: `rsquared` and `rsquared_adj` (the share of "
        "variation the model explains, with and without a penalty for the number "
        "of terms used) per (patient, treatment, feature) model, split by feature "
        "type and compartment. Most of the variation usually stays unexplained "
        "(`residual_pct = (1 - R²) × 100`), so expect low values."
    ),
    "UpSet (variate combinations)": (
        "A feature is a **hit** for a model term when `pvalue_fdr < threshold` "
        "and `coefficient > minimum` (increases only). We assign each feature "
        "to the set of terms it is a hit for; bars show how many features share "
        "each exact combination, which the viewer computes live. *Scope* sets the "
        "grouping (all models, per patient, per treatment, per tumor type, …). The "
        "overlap matrix below counts how many features every pair of groups shares."
    ),
    "Unique variates per model": (
        "A feature × (patient, treatment) yes/no heatmap of where treatment is a "
        "hit on its own versus together with at least one term from the chosen "
        "comparison group (biological or technical factors). Hit = `pvalue_fdr < "
        "threshold and coefficient > minimum`. The drill-down table below lists "
        "every term of one model and which ones are hits."
    ),
    # ---- 5.differential_analysis ----
    "Morphology vs viability": (
        "Does a bigger change in organoid appearance go along with a bigger change "
        "in viability? One point per patient, treatment, dose, and compartment. "
        "**x**: how different the treated organoids look from that patient's own "
        "DMSO control, averaged across features. **y**: how much viability changed "
        "relative to DMSO. We leave out points with no viability measurement."
    ),
    "Viability by tumor type": (
        "Does viability change differently by tumor type? We show viability as a "
        "percent of each patient's own DMSO control. The per-patient table has one "
        "row per patient and treatment; the summary table groups those by tumor "
        "type, with the number of patients and the mean, standard deviation, and "
        "median. MPNST and Other have only two patients each, so treat their "
        "averages with caution."
    ),
    "Plate position (inner vs outer)": (
        "Do wells on the outer ring of the plate (rows B and G, columns 2 and 11) "
        "differ from inner wells? Row B is only present on three of the twelve "
        "plates, so the row-B part of the outer ring comes from those three. "
        "We compare position only within plate x treatment "
        "x dose groups that contain both, because treatment and position overlap "
        "too much across the plate otherwise. Per feature, the effect is "
        "mean(outer) - mean(inner), averaged over groups. We compare the "
        "**global** effect (a root-mean-square summary across all features) "
        "against a chance baseline: we randomly reshuffle inner/outer well labels "
        "within each group 10,000 times to estimate that baseline. "
        "Per-feature p-values come from this same reshuffling approach; we then "
        "adjust them for testing many features at once (Benjamini-Hochberg q-values)."
    ),
    "DMSO column check": (
        "A sanity check for the plate-position result. DMSO occupies column 4 and "
        "column 9 of every plate (eight wells). If position does not matter, DMSO "
        "wells should be no more similar within a column than between columns. The "
        "statistic is the mean within-column correlation minus the mean "
        "between-column correlation, per patient plate; the comparison baseline comes "
        "from randomly reshuffling the column labels within each plate. The "
        "**all plates** row pools the plates."
    ),
    "Well correlation to DMSO": (
        "How closely each well matches DMSO: correlation between a well's profile "
        "and the median of the plate's *other* DMSO wells (we never compare a DMSO "
        "well to itself). Higher = more DMSO-like."
    ),
    "MEK feature tests": (
        "Per feature, this test compares the average MEK-inhibitor response (all "
        "four drugs, 1 and 10 µM pooled) against DMSO, Staurosporine, and Digoxin.\n"
        "- `mean`: average difference from the reference\n"
        "- `t` / `q`: significance of that difference across patients, FDR-adjusted\n"
        "- `rank`: features ordered by effect size (|t|)\n"
        "- `frac_patients_same_sign`: how many patients agree on direction"
    ),
    "MEK significance by level": (
        "How many features reach q < 0.05 at each level of the MEK tests. **Within "
        "patient**: over that patient's MEK conditions. **Within tumor type**: we pool "
        "patient x condition values, which is optimistic because conditions from one "
        "patient are not independent. **Across patients**: per-patient responses. "
        "**Across tumor types (shared)**: do the tumor-type means share the "
        "effect? **Across tumor types (differ)**: a one-way ANOVA (a test for "
        "comparing averages across more than two groups) asks whether the effect "
        "depends on tumor type. Below, the feature categories (Texture, "
        "Granularity, Intensity, ...) of the significant features."
    ),
}

SUBPANEL_BACKGROUND = {
    "corr_pairs": (
        "**Replicate/treatment-level (pairs).** Correlations between aggregate or "
        "consensus profiles (one row per replicate or treatment, summarized "
        "from single cells). We store it as a long table with one row per unique "
        "sample pair (upper triangle only, since correlation is symmetric), plus a "
        "companion table of patient/treatment/dose metadata for annotation. "
        "Files exist for 2D (3 slice strategies) and 3D (4 normalization variants)."
    ),
    "corr_per_patient": (
        "**Per-patient single-cell.** We compute correlation matrices of 3D "
        "single-cell feature-selected profiles separately for each patient (one "
        "row per cell before correlating, unlike the pairs view). Choose a "
        "normalization variant and patient. We include only variants that are "
        "truly single-cell."
    ),
    "pca_scree": (
        "**Scree plot.** Explained variance ratio of each principal component; "
        "we save it with the embeddings. A steep drop means a few components "
        "capture most of the variation."
    ),
    "lm_upset_profile": (
        "**Profile** sets the unit the model used: `organoid` (one row per "
        "organoid) or `sc` (one row per cell). **Model set**: `original` uses "
        "only treatment and count factors; `technical` adds location-related "
        "factors (plate position, cell x/y/z, depth in the image stack)."
    ),
    "lm_variates_drilldown": (
        "**Drill-down.** Pick one patient, treatment, and feature to see every "
        "term of that single model: its coefficient, adjusted p-value, R², "
        "and whether it passes the current hit thresholds."
    ),
    "plate_map": (
        "**Plate map.** One cell per well (rows A-H, columns 1-12). The screen uses "
        "rows B-G and columns 2-11, so the edge is empty. Row B is only present on "
        "three of the twelve plates, so most plate maps show rows C-G. Hover for the "
        "treatment."
    ),
    "mek_drug_consistency": (
        "**Consistency across MEK drugs.** The mean difference from DMSO for each MEK "
        "drug and dose, for the top features by |t| of the selected comparison. If "
        "the top features overlap, the columns should agree, rather than one "
        "compound driving the pattern."
    ),
    "mek_patient_values": (
        "**Per patient.** One point per patient for the chosen feature and "
        "comparison. A large mean that comes from one or two patients is weaker "
        "evidence than a consistent sign across patients."
    ),
}


def module_background(module: str) -> None:
    """Collapsed expander with the module-level background."""
    text = MODULE_BACKGROUND.get(module)
    if text:
        with st.expander(f"About {module}"):
            st.markdown(text)


def section_background(section: str) -> None:
    """Collapsed expander with the section-level background."""
    text = SECTION_BACKGROUND.get(section)
    if text:
        with st.expander(f"About: {section}"):
            st.markdown(text)


def subpanel_background(key: str, label: str = "About this panel") -> None:
    """Collapsed expander for a sub-panel inside a section."""
    text = SUBPANEL_BACKGROUND.get(key)
    if text:
        with st.expander(label):
            st.markdown(text)
