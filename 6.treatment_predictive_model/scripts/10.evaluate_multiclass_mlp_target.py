#!/usr/bin/env python
# coding: utf-8

# This notebook evaluates the global target/MOA multiclass MLP trained in `9.train_multiclass_mlp_target.ipynb`.
#
# It re-derives the held-out test-set predictions (same split, same deterministically-refit `StandardScaler`, saved model weights) and saves: per-cell predictions, overall and per-class accuracy/precision/recall/F1, and top-k accuracy. Plotting (including the confusion matrix) is done separately in `11.plot_multiclass_mlp_target.ipynb` (R).

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from notebook_init_utils import init_notebook
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    f1_score,
    top_k_accuracy_score,
)
from sklearn.preprocessing import LabelEncoder, StandardScaler
from torch.utils.data import DataLoader, TensorDataset

root_dir, in_notebook = init_notebook()

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
).resolve(strict=True)
model_results_dir = pathlib.Path(
    f"{root_dir}/6.treatment_predictive_model/model_results"
).resolve()
model_results_dir.mkdir(parents=True, exist_ok=True)


# ## Rebuild the same data/label setup used in training
# This mirrors `9.train_multiclass_mlp_target.ipynb` exactly, so refitting the `StandardScaler` and `LabelEncoder` on `train_df` reproduces what the saved model weights were trained under.

# In[3]:


LABEL_COL = "Metadata_Experiment_Target"

indices_df = pd.read_parquet(split_idx_file_path)
df = pd.read_parquet(profile_file_path)

feature_cols = [col for col in df.columns if not col.startswith("Metadata_")]

train_idx = indices_df["train_idx"].dropna().astype(int).tolist()
test_idx = indices_df["test_idx"].dropna().astype(int).tolist()
train_df, test_df = df.iloc[train_idx], df.iloc[test_idx]

label_encoder = LabelEncoder()
label_encoder.fit(train_df[LABEL_COL])
num_classes = len(label_encoder.classes_)

scaler = StandardScaler()
scaler.fit(train_df[feature_cols].values)
X_test = scaler.transform(test_df[feature_cols].values)
y_test = label_encoder.transform(test_df[LABEL_COL])

print(f"Test: {len(test_df)} rows, {num_classes} classes")


# ## Model architecture
# Must match `9.train_multiclass_mlp_target.ipynb` exactly so the saved `state_dict` loads correctly.

# In[4]:


HIDDEN_DIMS = (128, 64)


class MulticlassMLP(nn.Module):
    """A feed-forward multiclass classifier over single-cell morphology features."""

    def __init__(
        self, input_dim: int, num_classes: int, hidden_dims: tuple = HIDDEN_DIMS
    ):
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
        layers.append(nn.Linear(prev_dim, num_classes))
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x)


model = MulticlassMLP(input_dim=len(feature_cols), num_classes=num_classes).to(device)
model.load_state_dict(
    torch.load(model_output_dir / "multiclass_mlp_target.pt", map_location=device)
)
model.eval()


# ## Score the test set

# In[5]:


test_loader = DataLoader(
    TensorDataset(
        torch.tensor(X_test, dtype=torch.float32),
        torch.tensor(y_test, dtype=torch.long),
    ),
    batch_size=128,
    shuffle=False,
)

all_logits, all_labels = [], []
with torch.no_grad():
    for X_batch, y_batch in test_loader:
        all_logits.append(model(X_batch.to(device)).cpu())
        all_labels.append(y_batch)
logits = torch.cat(all_logits).numpy()
labels = torch.cat(all_labels).numpy()
preds = logits.argmax(axis=1)
probs = torch.softmax(torch.tensor(logits), dim=1).numpy()

test_predictions_df = pd.DataFrame(
    {
        "Metadata_true_label": label_encoder.inverse_transform(labels),
        "Metadata_predicted_label": label_encoder.inverse_transform(preds),
    }
)
test_predictions_df.to_parquet(
    model_results_dir / "multiclass_mlp_target_test_predictions.parquet",
    index=False,
)
test_predictions_df.head()


# ## Overall metrics and top-k accuracy

# In[6]:


top_k_values = [1, 2, 3, 5]
top_k_rows = [
    {
        "k": k,
        "top_k_accuracy": top_k_accuracy_score(
            labels, probs, k=k, labels=np.arange(num_classes)
        ),
    }
    for k in top_k_values
]
top_k_df = pd.DataFrame(top_k_rows)
top_k_df.to_parquet(
    model_results_dir / "multiclass_mlp_target_top_k_accuracy.parquet", index=False
)

print(f"Accuracy: {accuracy_score(labels, preds):.4f}")
print(f"Macro F1: {f1_score(labels, preds, average='macro'):.4f}")
print(f"Weighted F1: {f1_score(labels, preds, average='weighted'):.4f}")
top_k_df


# ## Per-class precision/recall/F1

# In[7]:


report_dict = classification_report(
    labels,
    preds,
    target_names=label_encoder.classes_,
    output_dict=True,
    zero_division=0,
)
report_df = (
    pd.DataFrame(report_dict)
    .T.reset_index(names="Metadata_class")
    .query("Metadata_class not in ['accuracy']")
)
report_df.to_parquet(
    model_results_dir / "multiclass_mlp_target_classification_report.parquet",
    index=False,
)
report_df.sort_values("f1-score", ascending=False)
