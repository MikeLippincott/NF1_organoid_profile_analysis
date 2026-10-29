#!/usr/bin/env python
# coding: utf-8

# The goal of this notebook is to train a per-patient model that classifies cells as DMSO-like (negative control) or cytotoxic-like (positive control) based on their image-based morphological profiles.
# Applying this model to treated cells lets us measure how similar each cell is to the negative and positive controls.
# We are particularly interested in cells that differ from DMSO but do not resemble the cytotoxic phenotype, as these may reflect treatment-specific, non-toxic responses.
# To identify these cells, we [describe approach, e.g., combine the classifier probability with each cell's distance from the DMSO profile].

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from notebook_init_utils import init_notebook
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from tqdm.notebook import tqdm as tqdm_notebook

root_dir, in_notebook = init_notebook()
progress_bar = tqdm_notebook if in_notebook else tqdm

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.manual_seed(0)


# In[2]:


profile_file_path = pathlib.Path(
    f"{root_dir}/data/profiles_3D/all_patients/1.feature_selected_profiles/sc_norm_fs_profiles.parquet"
).resolve(strict=True)

split_idx_file_path = pathlib.Path(
    f"{root_dir}/6.treatment_predictive_model/data_splits/sc_norm_fs_profiles_split_idx.parquet"
).resolve(strict=True)

model_output_dir = pathlib.Path(
    f"{root_dir}/6.treatment_predictive_model/trained_models"
).resolve()
model_results_dir = pathlib.Path(
    f"{root_dir}/6.treatment_predictive_model/model_results"
).resolve()
model_output_dir.mkdir(parents=True, exist_ok=True)
model_results_dir.mkdir(parents=True, exist_ok=True)


# In[3]:


indices_df = pd.read_parquet(split_idx_file_path)
df = pd.read_parquet(profile_file_path)


# In[4]:


# get the train, validation, and test indices from the dataframe
train_idx = indices_df["train_idx"].dropna().astype(int).tolist()
val_idx = indices_df["val_idx"].dropna().astype(int).tolist()
test_idx = indices_df["test_idx"].dropna().astype(int).tolist()

train_df, val_df, test_df = df.iloc[train_idx], df.iloc[val_idx], df.iloc[test_idx]


# ## Define the two classes
# `DMSO` is the negative control (vehicle only, "normal" morphology) and `Staurosporine` is the cytotoxic positive control. Every other treatment is excluded from training - the resulting classifier will later be applied to those treated cells to score how DMSO-like vs. cytotoxic-like they appear.
#
# A separate MLP is trained per patient, since patients are expected to differ in their baseline (DMSO) and cytotoxic-response morphology.

# In[5]:


PATIENT_COL = "Metadata_Biology_PatientTumor"
TREATMENT_COL = "Metadata_Experiment_Treatment"
LABEL_COL = "Metadata_control_label"
CONTROL_LABEL_MAP = {"DMSO": 0, "Staurosporine": 1}

feature_cols = [col for col in df.columns if not col.startswith("Metadata_")]

train_df = train_df[train_df[TREATMENT_COL].isin(CONTROL_LABEL_MAP)].copy()
val_df = val_df[val_df[TREATMENT_COL].isin(CONTROL_LABEL_MAP)].copy()
test_df = test_df[test_df[TREATMENT_COL].isin(CONTROL_LABEL_MAP)].copy()

for split_df in (train_df, val_df, test_df):
    split_df[LABEL_COL] = (
        split_df[TREATMENT_COL].map(CONTROL_LABEL_MAP).astype(np.float32)
    )

print(
    f"Train: {len(train_df)} rows, Val: {len(val_df)} rows, Test: {len(test_df)} rows"
)
print(f"Number of features: {len(feature_cols)}")


# ## MLP model and training/evaluation helpers

# In[6]:


BATCH_SIZE = 64
MAX_EPOCHS = 100
PATIENCE = 10
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-4
HIDDEN_DIMS = (128, 64)


class MLPClassifier(nn.Module):
    """A feed-forward binary classifier over single-cell morphology features."""

    def __init__(self, input_dim: int, hidden_dims: tuple = HIDDEN_DIMS):
        super().__init__()
        layers = []
        prev_dim = input_dim
        for hidden_dim in hidden_dims:
            layers += [
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.3),
            ]
            prev_dim = hidden_dim
        layers.append(nn.Linear(prev_dim, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).squeeze(-1)


def make_loader(
    X: np.ndarray,
    y: np.ndarray,
    batch_size: int,
    shuffle: bool,
    drop_last: bool = False,
) -> DataLoader:
    dataset = TensorDataset(
        torch.tensor(X, dtype=torch.float32), torch.tensor(y, dtype=torch.float32)
    )
    # drop_last avoids handing BatchNorm1d a final batch of size 1 while
    # training (it requires >1 sample per channel to compute batch statistics).
    return DataLoader(
        dataset, batch_size=batch_size, shuffle=shuffle, drop_last=drop_last
    )


def train_mlp(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader,
    pos_weight: torch.Tensor | None,
) -> tuple:
    optimizer = torch.optim.Adam(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    best_val_loss = float("inf")
    best_state = None
    epochs_without_improvement = 0
    history = {"train_loss": [], "val_loss": []}

    for _ in range(MAX_EPOCHS):
        model.train()
        running_loss = 0.0
        n_train_seen = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            loss = criterion(model(X_batch), y_batch)
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * X_batch.size(0)
            n_train_seen += X_batch.size(0)
        train_loss = running_loss / n_train_seen

        model.eval()
        val_running_loss = 0.0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                val_running_loss += criterion(
                    model(X_batch), y_batch
                ).item() * X_batch.size(0)
        val_loss = val_running_loss / len(val_loader.dataset)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= PATIENCE:
                break

    model.load_state_dict(best_state)
    return model, history


def evaluate_mlp(model: nn.Module, data_loader: DataLoader) -> dict:
    model.eval()
    all_logits, all_labels = [], []
    with torch.no_grad():
        for X_batch, y_batch in data_loader:
            all_logits.append(model(X_batch.to(device)).cpu())
            all_labels.append(y_batch)
    logits = torch.cat(all_logits).numpy()
    labels = torch.cat(all_labels).numpy()
    probs = 1 / (1 + np.exp(-logits))
    preds = (probs >= 0.5).astype(int)

    metrics = {
        "accuracy": accuracy_score(labels, preds),
        "f1": f1_score(labels, preds, zero_division=0),
        "roc_auc": (
            roc_auc_score(labels, probs) if len(np.unique(labels)) > 1 else np.nan
        ),
    }
    return metrics


# ## Train one MLP per patient
# Features are standardized per patient using only that patient's training rows. Patients without both classes present in every split (e.g. no `Staurosporine` wells) are skipped.

# In[7]:


patients = sorted(train_df[PATIENT_COL].unique())
per_patient_metrics = []
per_patient_histories = {}

for patient in progress_bar(patients, desc="Training per-patient MLPs"):
    patient_train = train_df[train_df[PATIENT_COL] == patient]
    patient_val = val_df[val_df[PATIENT_COL] == patient]
    patient_test = test_df[test_df[PATIENT_COL] == patient]

    if (
        patient_train[LABEL_COL].nunique() < 2
        or patient_val[LABEL_COL].nunique() < 2
        or patient_test[LABEL_COL].nunique() < 2
    ):
        print(f"Skipping {patient}: missing a class in train/val/test.")
        continue

    scaler = StandardScaler()
    X_train = scaler.fit_transform(patient_train[feature_cols].values)
    X_val = scaler.transform(patient_val[feature_cols].values)
    X_test = scaler.transform(patient_test[feature_cols].values)

    y_train = patient_train[LABEL_COL].values
    y_val = patient_val[LABEL_COL].values
    y_test = patient_test[LABEL_COL].values

    # drop_last=True only for training, and only when it leaves at least one
    # batch - otherwise BatchNorm1d can be handed a final training batch of
    # size 1, which it cannot compute batch statistics from.
    train_drop_last = len(X_train) > BATCH_SIZE
    train_loader = make_loader(
        X_train, y_train, BATCH_SIZE, shuffle=True, drop_last=train_drop_last
    )
    val_loader = make_loader(X_val, y_val, BATCH_SIZE, shuffle=False)
    test_loader = make_loader(X_test, y_test, BATCH_SIZE, shuffle=False)

    n_pos = y_train.sum()
    n_neg = len(y_train) - n_pos
    pos_weight = torch.tensor([n_neg / n_pos], dtype=torch.float32, device=device)

    model = MLPClassifier(input_dim=len(feature_cols)).to(device)
    model, history = train_mlp(model, train_loader, val_loader, pos_weight)
    per_patient_histories[patient] = history

    test_metrics = evaluate_mlp(model, test_loader)
    test_metrics[PATIENT_COL] = patient
    test_metrics["n_train"] = len(patient_train)
    test_metrics["n_val"] = len(patient_val)
    test_metrics["n_test"] = len(patient_test)
    test_metrics["n_epochs_trained"] = len(history["train_loss"])
    per_patient_metrics.append(test_metrics)

    model_output_path = model_output_dir / f"binary_mlp_control_{patient}.pt"
    torch.save(model.state_dict(), model_output_path)

metrics_df = pd.DataFrame(per_patient_metrics)
metrics_df


# In[8]:


metrics_output_path = model_results_dir / "binary_mlp_control_test_metrics.parquet"
metrics_df.to_parquet(metrics_output_path, index=False)


# ## Save per-patient training curves
# Plotting (per-patient test metrics, training curves) is done separately in `3.plot_binary_mlp_control_performance.ipynb` (R) - this notebook only computes and saves the underlying data.

# In[9]:


history_records = []
for patient, history in per_patient_histories.items():
    history_records.append(
        pd.DataFrame(
            {
                PATIENT_COL: patient,
                "epoch": range(1, len(history["train_loss"]) + 1),
                "train_loss": history["train_loss"],
                "val_loss": history["val_loss"],
            }
        )
    )

history_df = pd.concat(history_records, ignore_index=True)
history_output_path = model_results_dir / "binary_mlp_control_training_curves.parquet"
history_df.to_parquet(history_output_path, index=False)
history_df.head()
