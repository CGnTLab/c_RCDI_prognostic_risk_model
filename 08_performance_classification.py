#!/usr/bin/env python
# coding: utf-8


# Install required packages (if not already installed)
#!pip install imbalanced-learn xgboost scikit-learn seaborn matplotlib --quiet

#!pip install -U imbalanced-learn

# Imports
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, RandomizedSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, auc, confusion_matrix
from sklearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline as ImbPipeline
import matplotlib.pyplot as plt
import seaborn as sns
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (classification_report, roc_auc_score, roc_curve, auc,
                             confusion_matrix, accuracy_score, precision_score,
                             recall_score, f1_score, matthews_corrcoef)
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.metrics import make_scorer

# load datasets
train = pd.read_csv("../Train_RCDI_categorized.csv", index_col = 0)
test = pd.read_csv("../Test_RCDI_categorized.csv", index_col = 0)
val1 = pd.read_csv("../GSE39582_RCDI_categorized.csv", index_col = 0)
val2 = pd.read_csv("../GSE161158_RCDI_categorized.csv", index_col = 0)
train.drop(columns = ["OS_MONTHS","OS_STATUS","RCDI"], inplace =True)
test.drop(columns = ["OS_MONTHS","OS_STATUS","RCDI"], inplace=True)
val1.drop(columns = ["OS_MONTHS","OS_STATUS","RCDI"], inplace=True)
val2.drop(columns = ["OS_MONTHS","OS_STATUS","RCDI"], inplace=True)

train.columns, test.columns, val1.columns, val2.columns


# Features and Labels
X_train = train.drop(columns=["RiskClass"])
y_train = train["RiskClass"].map({"low": 0, "high": 1})

X_test = test.drop(columns=["RiskClass"])
y_test = test["RiskClass"].map({"low": 0, "high": 1})

X_val1 = val1.drop(columns=["RiskClass"])
y_val1 = val1["RiskClass"].map({"low": 0, "high": 1})

X_val2 = val2.drop(columns=["RiskClass"])
y_val2 = val2["RiskClass"].map({"low": 0, "high": 1})

X_train.shape, y_train.shape, X_test.shape, y_test.shape, X_val1.shape, y_val1.shape, X_val2.shape, y_val2.shape


# Standardization
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)
X_val1 = scaler.transform(X_val1)
X_val2 = scaler.transform(X_val2)


# Define classifiers and hyperparameter grids
classifiers = {
    "LogisticRegression": (LogisticRegression(max_iter=1000, random_state=42), {
        "clf__C": np.logspace(-3, 3, 10),
        "clf__penalty": ["l2"]
    }),
    "SVM": (SVC(probability=True, random_state=42), {
        "clf__C": [0.1, 1, 10, 100],
        "clf__kernel": ["linear", "rbf"]
    }),
    "RandomForest": (RandomForestClassifier(random_state=42), {
        "clf__n_estimators": [100, 200, 500],
        "clf__max_depth": [None, 10, 20],
        "clf__min_samples_split": [2, 5, 10]
    }),
    "XGBoost": (XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42), {
        "clf__n_estimators": [100, 200],
        "clf__max_depth": [3, 5, 7],
        "clf__learning_rate": [0.01, 0.1, 0.2]
    }),
    "KNN": (KNeighborsClassifier(), {
        "clf__n_neighbors": [3, 5, 7, 9]
    }),
    "NaiveBayes": (GaussianNB(), {})
}

## Use cross validation for training dataset
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.base import ClassifierMixin

# containers
results = {}
roc_data = {}
best_models = {}
val_results = {}
train_roc_data = {}
performance_table = []

cv10 = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
        
# Define scoring metrics
scoring = {
    'accuracy': 'accuracy',
    'precision': 'precision',
    'recall': 'recall',
    'f1': 'f1',
    'roc_auc': 'roc_auc',
    'mcc': make_scorer(matthews_corrcoef)
}

for name, (clf, params) in classifiers.items():
    print(f"\nTraining {name}...")

    try:
        if name == "XGBoost":
            # Apply SMOTE externally
            sm = SMOTE(random_state=42)
            X_train_res, y_train_res = sm.fit_resample(X_train, y_train)

            pipe = ImbPipeline([("clf", clf)])
            search = RandomizedSearchCV(pipe, param_distributions=params, scoring="roc_auc", cv=cv10, n_iter=10, n_jobs=-1, random_state=42)
            search.fit(X_train_res, y_train_res)
            best_model = search.best_estimator_

            # Manually compute CV metrics
            accs, precs, recs, f1s, mccs, aucs = [], [], [], [], [], []

            for train_idx, val_idx in cv10.split(X_train_res, y_train_res):
                X_tr, X_val = X_train_res[train_idx], X_train_res[val_idx]
                y_tr, y_val = y_train_res[train_idx], y_train_res[val_idx]

                best_model.fit(X_tr, y_tr)
                y_val_pred = best_model.predict(X_val)
                y_val_proba = best_model.predict_proba(X_val)[:, 1]

                accs.append(accuracy_score(y_val, y_val_pred))
                precs.append(precision_score(y_val, y_val_pred, zero_division=0))
                recs.append(recall_score(y_val, y_val_pred, zero_division=0))
                f1s.append(f1_score(y_val, y_val_pred, zero_division=0))
                mccs.append(matthews_corrcoef(y_val, y_val_pred))
                if len(np.unique(y_val)) == 2:
                    aucs.append(roc_auc_score(y_val, y_val_proba))

            cv_results = {
                'test_accuracy': accs,
                'test_precision': precs,
                'test_recall': recs,
                'test_f1': f1s,
                'test_mcc': mccs,
                'test_roc_auc': aucs
            }

            train_auc_mean = np.mean(aucs) if aucs else np.nan
            train_auc_std = np.std(aucs) if aucs else np.nan

        else:
            pipe = ImbPipeline([
                ("smote", SMOTE(random_state=42)),
                ("clf", clf)
            ])
            search = RandomizedSearchCV(pipe, param_distributions=params, scoring="roc_auc", cv=cv10, n_iter=10, n_jobs=-1, random_state=42)
            search.fit(X_train, y_train)
            best_model = search.best_estimator_

            # Standard cross-validation
            cv_results = cross_validate(
                best_model, X_train, y_train,
                scoring=scoring, cv=cv10, n_jobs=-1
            )

            auc_values = cv_results.get('test_roc_auc', [])
            train_auc_mean = np.mean(auc_values) if len(auc_values) > 0 else np.nan
            train_auc_std = np.std(auc_values) if len(auc_values) > 0 else np.nan

        best_models[name] = best_model
        train_roc_data[name] = train_auc_mean

        print(f"{name} - 10-Fold Cross-Validated AUC: Mean = {train_auc_mean:.4f}, Std = {train_auc_std:.4f}")

        # Test set evaluation
        y_test_pred = best_model.predict(X_test)
        y_test_proba = best_model.predict_proba(X_test)[:, 1]
        test_auc = roc_auc_score(y_test, y_test_proba)
        fpr, tpr, _ = roc_curve(y_test, y_test_proba)
        roc_data[name] = (fpr, tpr)

        results[name] = {
            "Train AUC": train_auc_mean,
            "Test AUC": test_auc
        }

        # Validation sets
        val_metrics = {}
        for val_name, (X_val, y_val) in zip(["GSE39582", "GSE161158"], [(X_val1, y_val1), (X_val2, y_val2)]):
            y_val_pred = best_model.predict(X_val)
            y_val_proba = best_model.predict_proba(X_val)[:, 1]
            val_auc = roc_auc_score(y_val, y_val_proba)
            val_results[(name, val_name)] = val_auc
            val_metrics[val_name] = {
                "accuracy": accuracy_score(y_val, y_val_pred),
                "precision": precision_score(y_val, y_val_pred),
                "recall": recall_score(y_val, y_val_pred),
                "f1": f1_score(y_val, y_val_pred),
                "mcc": matthews_corrcoef(y_val, y_val_pred),
                "auc": val_auc
            }

        # Add to performance table
        performance_table.append({
            "Model": name,

            # Train (from CV)
            "Train Accuracy": np.mean(cv_results['test_accuracy']),
            "Train Precision": np.mean(cv_results['test_precision']),
            "Train Recall": np.mean(cv_results['test_recall']),
            "Train F1": np.mean(cv_results['test_f1']),
            "Train MCC": np.mean(cv_results['test_mcc']),
            "Train AUC": train_auc_mean,

            # Test
            "Test Accuracy": accuracy_score(y_test, y_test_pred),
            "Test Precision": precision_score(y_test, y_test_pred),
            "Test Recall": recall_score(y_test, y_test_pred),
            "Test F1": f1_score(y_test, y_test_pred),
            "Test MCC": matthews_corrcoef(y_test, y_test_pred),
            "Test AUC": test_auc,

            # Validation sets
            "GSE39582 Accuracy": val_metrics["GSE39582"]["accuracy"],
            "GSE39582 Precision": val_metrics["GSE39582"]["precision"],
            "GSE39582 Recall": val_metrics["GSE39582"]["recall"],
            "GSE39582 F1": val_metrics["GSE39582"]["f1"],
            "GSE39582 MCC": val_metrics["GSE39582"]["mcc"],
            "GSE39582 AUC": val_metrics["GSE39582"]["auc"],

            "GSE161158 Accuracy": val_metrics["GSE161158"]["accuracy"],
            "GSE161158 Precision": val_metrics["GSE161158"]["precision"],
            "GSE161158 Recall": val_metrics["GSE161158"]["recall"],
            "GSE161158 F1": val_metrics["GSE161158"]["f1"],
            "GSE161158 MCC": val_metrics["GSE161158"]["mcc"],
            "GSE161158 AUC": val_metrics["GSE161158"]["auc"],

            # AUC stats and best parameters
            "CV10 AUC Mean": train_auc_mean,
            "CV10 AUC Std": train_auc_std,
            "Best Params": search.best_params_
        })

    except Exception as e:
        print(f"❌ Error in model {name}: {str(e)}")
        train_roc_data[name] = np.nan


pd.DataFrame(performance_table)

df_perf = pd.DataFrame(performance_table)


# Convert from wide to long format for easier plotting
plot_df = pd.melt(
    df_perf,
    id_vars=["Model"],
    value_vars=[
        "Train Accuracy", "Test Accuracy", "GSE39582 Accuracy", "GSE161158 Accuracy",
        "Train Recall", "Test Recall", "GSE39582 Recall", "GSE161158 Recall"
    ],
    var_name="Metric_Dataset",
    value_name="Score"
)

# Split Metric_Dataset into two columns
plot_df[['Dataset', 'Metric']] = plot_df['Metric_Dataset'].str.extract(r'(\w+)\s+(\w+)')
plot_df.drop(columns='Metric_Dataset', inplace=True)


## Plot accuracy - recall barplots
import seaborn as sns
import matplotlib.pyplot as plt

# Set style
sns.set(style="whitegrid", context="talk")

for metric in ['Accuracy', 'Recall']:
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=plot_df[plot_df['Metric'] == metric],
        x="Model", y="Score", hue="Dataset",
        palette="Set2"
    )
    plt.title(f"{metric} Score Across Datasets for Each Classifier", fontsize=16)
    plt.ylabel(f"{metric} Score", fontsize=14)
    plt.xlabel("Classifier", fontsize=14)
    plt.ylim(0, 1)
    plt.legend(title="Dataset", loc='best', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

## Save the performance table
performance_table = pd.DataFrame(performance_table)

performance_table.to_csv("model_performances.csv", index=False)
performance_table.head()


## Create AURCO plots for the four datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import numpy as np
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE

cv10 = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
roc_cv_train = {}
mean_fpr = np.linspace(0, 1, 100)

for name, model in best_models.items():
    print(f"CV ROC for {name}...")
    tprs = []
    aucs = []

    if name == "XGBoost":
        X_train_res, y_train_res = SMOTE(random_state=42).fit_resample(X_train, y_train)
        X_data, y_data = X_train_res, y_train_res
    else:
        X_data, y_data = X_train, y_train

    for train_idx, val_idx in cv10.split(X_data, y_data):
        X_tr, X_val = X_data[train_idx], X_data[val_idx]
        y_tr, y_val = y_data[train_idx], y_data[val_idx]

        model.fit(X_tr, y_tr)
        y_proba = model.predict_proba(X_val)[:, 1]
        fpr, tpr, _ = roc_curve(y_val, y_proba)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(auc(fpr, tpr))

    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)

    roc_cv_train[name] = (mean_fpr, mean_tpr, std_tpr, mean_auc, std_auc)

def plot_mean_roc_curves_cv(
    roc_dict, title, filename, figsize=(6.5, 6), fontsize=13
):
    plt.figure(figsize=figsize)
    for model_name, (fpr, mean_tpr, _, mean_auc, _) in roc_dict.items():
        plt.plot(fpr, mean_tpr, lw=2,
                 label=f"{model_name} (AUC = {mean_auc:.2f})")

    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate", fontsize=fontsize)
    plt.ylabel("True Positive Rate", fontsize=fontsize)
    plt.title(title, fontsize=fontsize + 2)
    plt.legend(loc="lower right", fontsize=fontsize - 1)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(filename, dpi=600)
    plt.close()


def compute_standard_roc(models_dict, X, y):
    roc_dict = {}
    for name, model in models_dict.items():
        y_proba = model.predict_proba(X)[:, 1]
        fpr, tpr, _ = roc_curve(y, y_proba)
        auc_score = auc(fpr, tpr)
        roc_dict[name] = (fpr, tpr, auc_score)
    return roc_dict

def plot_standard_roc(roc_dict, title, filename, figsize=(6.5, 6), fontsize=13):
    plt.figure(figsize=figsize)
    for model_name, (fpr, tpr, auc_score) in roc_dict.items():
        plt.plot(fpr, tpr, lw=2, label=f"{model_name} (AUC = {auc_score:.2f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate", fontsize=fontsize)
    plt.ylabel("True Positive Rate", fontsize=fontsize)
    plt.title(title, fontsize=fontsize + 2)
    plt.legend(loc="lower right", fontsize=fontsize - 1)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(filename, dpi=600)
    plt.close()


# Cross-validated (Train)
plot_mean_roc_curves_cv(
    roc_cv_train,
    "Cross-Validated ROC Curve (Train Set)",
    "roc_train_cv.png"
)

# Standard Test Set
roc_test = compute_standard_roc(best_models, X_test, y_test)
plot_standard_roc(roc_test, "ROC Curve - Test Set", "roc_test.png")

# Validation Set 1 (GSE39582)
roc_val1 = compute_standard_roc(best_models, X_val1, y_val1)
plot_standard_roc(roc_val1, "ROC Curve - GSE39582", "roc_gse39582.png")

# Validation Set 2 (GSE161158)
roc_val2 = compute_standard_roc(best_models, X_val2, y_val2)
plot_standard_roc(roc_val2, "ROC Curve - GSE161158", "roc_gse161158.png")



# Create summary DataFrame for AUC scores
summary_auc = pd.DataFrame(columns=["Train", "Test", "GSE39582", "GSE161158"])

for model_name in best_models.keys():
    # Get CV-based training AUC (mean from cross_validate)
    train_auc = train_roc_data[model_name]
    
    # Get test AUC (directly computed)
    test_auc = results[model_name]["Test AUC"]
    
    # Get validation AUCs
    gse39582_auc = val_results.get((model_name, "GSE39582"), None)
    gse161158_auc = val_results.get((model_name, "GSE161158"), None)
    
    summary_auc.loc[model_name] = [train_auc, test_auc, gse39582_auc, gse161158_auc]

# Round for clarity
summary_auc = summary_auc.round(3)

# Save to CSV
summary_auc.to_csv("classifier_auc_summary.csv")

# Display
print(summary_auc)


# Feature importance plots for RF and XGB
for name in ["RandomForest", "XGBoost"]:
    model = best_models[name].named_steps['clf']
    importances = model.feature_importances_
    feat_names = train.drop(columns=["RiskClass"]).columns
    importance_df = pd.DataFrame({"Feature": feat_names, "Importance": importances}).sort_values(by="Importance", ascending=False).head(20)
    plt.figure(figsize=(14, 10))
    sns.barplot(x="Importance", y="Feature", data=importance_df)
    plt.title(f"All Feature Importances ({name})")
    plt.tight_layout()
    plt.savefig(f"feature_importance_{name}.png", dpi=600)
    plt.show()

# CONSOLIDATED PERMUTATION AND SHAP FEATURE IMPORTANCE ANALYSIS

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.inspection import permutation_importance
import shap
from imblearn.over_sampling import SMOTE
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.titlesize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.format': 'svg',
    'figure.figsize': (3.5, 2.8),
    'axes.linewidth': 0.5,
    'lines.linewidth': 1.0,
    'patch.linewidth': 0.5
})

feature_names = train.drop(columns=["RiskClass"]).columns.tolist()

def normalize_per_model(df, model_col="Model", imp_col="Importance"):
    df = df.copy()
    for model in df[model_col].unique():
        mask = df[model_col] == model
        vals = df.loc[mask, imp_col]
        minv, maxv = vals.min(), vals.max()
        if maxv > minv:
            df.loc[mask, imp_col] = (vals - minv) / (maxv - minv)
    return df

def compute_permutation_all(models_dict, X, y, feature_names, scoring='roc_auc', n_repeats=10, sample_size=1000):
    from sklearn.base import BaseEstimator, ClassifierMixin

    X_eval, _, y_eval, _ = train_test_split(
        X, y,
        train_size=min(sample_size, len(X)),
        stratify=y,
        random_state=42
    )
    all_results = []

    class XGBWrapper(BaseEstimator, ClassifierMixin):
        """Robust sklearn-compatible wrapper for XGBoost classifiers."""
        def __init__(self, xgb_model):
            # Ensure we store the actual model, not a dict or other container
            if hasattr(xgb_model, 'predict') and hasattr(xgb_model, 'fit'):
                self.xgb_model = xgb_model
            else:
                raise ValueError(f"Expected XGBoost model with predict/fit methods, got {type(xgb_model)}")

        def fit(self, X, y, **kwargs):
            self.xgb_model.fit(X, y, **kwargs)
            return self

        def predict(self, X):
            return self.xgb_model.predict(X)

        def predict_proba(self, X):
            if hasattr(self.xgb_model, "predict_proba"):
                return self.xgb_model.predict_proba(X)
            else:
                # Fallback for models that don't have predict_proba
                preds = self.xgb_model.predict(X)
                if preds.ndim == 1:
                    return np.column_stack([1 - preds, preds])
                return preds

        def get_params(self, deep=True):
            if hasattr(self.xgb_model, 'get_params'):
                return self.xgb_model.get_params(deep=deep)
            return {}

        def set_params(self, **params):
            if hasattr(self.xgb_model, 'set_params'):
                self.xgb_model.set_params(**params)
            return self

        @property
        def classes_(self):
            if hasattr(self.xgb_model, 'classes_'):
                return self.xgb_model.classes_
            return np.array([0, 1])  # Default binary classes

        def _estimator_type(self):
            return "classifier"

        def __sklearn_tags__(self):
            return {
                "binary_only": False,
                "multilabel": False,
                "multioutput": False,
                "requires_y": True,
                "classifier": True,
                "regressor": False,
            }

    def extract_final_estimator(model):
        """Recursively extract the trained estimator from pipelines or nested structures."""
        if model is None:
            return None
            
        # Handle pipeline structures
        if hasattr(model, "named_steps"):
            if "clf" in model.named_steps:
                return extract_final_estimator(model.named_steps["clf"])
            # Try other common step names
            for step_name in ["classifier", "estimator", "model"]:
                if step_name in model.named_steps:
                    return extract_final_estimator(model.named_steps[step_name])
                    
        # Handle GridSearchCV/RandomizedSearchCV
        if hasattr(model, "best_estimator_"):
            return extract_final_estimator(model.best_estimator_)
            
        # Handle dictionary storage (your case)
        if isinstance(model, dict):
            # Look for model-like objects in the dictionary
            for key, value in model.items():
                if hasattr(value, "predict") and hasattr(value, "fit"):
                    return extract_final_estimator(value)
                elif hasattr(value, "named_steps"):
                    return extract_final_estimator(value)
                elif hasattr(value, "best_estimator_"):
                    return extract_final_estimator(value)
                    
        # If it has predict method, assume it's the final estimator
        if hasattr(model, "predict") and hasattr(model, "fit"):
            return model
            
        return None

    for name, model in models_dict.items():
        print(f"Calculating permutation importance for {name}...")
        try:
            clf = extract_final_estimator(model)
            if clf is None:
                print(f"  Could not extract valid estimator from {name}")
                continue

            # Debug: print what we extracted
            print(f"  Extracted estimator type: {type(clf)}")

            # For XGBoost, wrap it properly
            model_name_lower = name.lower()
            if "xgb" in model_name_lower or "xgboost" in model_name_lower:
                print(f"  Wrapping {name} with XGBWrapper")
                try:
                    clf = XGBWrapper(clf)
                    print(f"  Successfully wrapped, wrapper type: {type(clf)}")
                except Exception as wrap_ex:
                    print(f"  Failed to wrap XGBoost: {wrap_ex}")
                    continue

            # Ensure clf is not a dict before proceeding
            if isinstance(clf, dict):
                print(f"  Error: Extracted object is still a dict for {name}")
                continue

            result = permutation_importance(
                clf,
                X_eval,
                y_eval,
                scoring=scoring,
                n_repeats=n_repeats,
                random_state=42,
                n_jobs=-1
            )

            df = pd.DataFrame({
                'Gene': feature_names,
                'Importance': result.importances_mean,
                'Std': result.importances_std,
                'Model': name
            })
            all_results.append(df)
            print(f"  ✓ Completed {name}")

        except Exception as ex:
            print(f"  ✗ Error with {name}: {ex}")
            import traceback
            traceback.print_exc()

    if not all_results:
        raise RuntimeError("No permutation results computed.")
    result = pd.concat(all_results, axis=0)
    return normalize_per_model(result)

def compute_shap_all(models_dict, X, feature_names, sample_size=200, background_size=20):
    X_sample = X if len(X) <= sample_size else X[np.random.choice(len(X), sample_size, replace=False)]
    all_results = []
    
    def extract_final_estimator(model):
        """Same extraction logic as above"""
        if model is None:
            return None
        if hasattr(model, "named_steps"):
            if "clf" in model.named_steps:
                return extract_final_estimator(model.named_steps["clf"])
            for step_name in ["classifier", "estimator", "model"]:
                if step_name in model.named_steps:
                    return extract_final_estimator(model.named_steps[step_name])
        if hasattr(model, "best_estimator_"):
            return extract_final_estimator(model.best_estimator_)
        if isinstance(model, dict):
            for key, value in model.items():
                if hasattr(value, "predict") and hasattr(value, "fit"):
                    return extract_final_estimator(value)
                elif hasattr(value, "named_steps"):
                    return extract_final_estimator(value)
                elif hasattr(value, "best_estimator_"):
                    return extract_final_estimator(value)
        if hasattr(model, "predict") and hasattr(model, "fit"):
            return model
        return None
    
    for name, model in models_dict.items():
        print(f"Calculating SHAP importance for {name}...")
        try:
            clf = extract_final_estimator(model)
            if clf is None:
                print(f"  Could not extract valid estimator from {name}")
                continue
                
            print(f"  Extracted estimator type: {type(clf)}")
            
            # Choose best explainer per model
            if "RandomForest" in name or "XGBoost" in name:
                explainer = shap.TreeExplainer(clf)
                shap_values = explainer.shap_values(X_sample)
            elif "LogisticRegression" in name:
                background = shap.sample(X, min(background_size, len(X)))
                explainer = shap.LinearExplainer(clf, background, feature_perturbation="interventional")
                shap_values = explainer.shap_values(X_sample)
            else:
                # fallback for SVM, KNN, NaiveBayes
                if len(X) > background_size:
                    kmeans = KMeans(n_clusters=background_size, random_state=42)
                    kmeans.fit(X)
                    background = kmeans.cluster_centers_
                else:
                    background = X
                explainer = shap.KernelExplainer(clf.predict_proba, background, approximate=True)
                shap_values = explainer.shap_values(X_sample)

            # Robust aggregation
            if isinstance(shap_values, list):
                shap_abs = np.mean([np.abs(sv) for sv in shap_values], axis=0)
            elif shap_values.ndim == 3:
                shap_abs = np.mean(np.abs(shap_values), axis=(0, 2))
            elif shap_values.ndim == 2:
                shap_abs = np.mean(np.abs(shap_values), axis=0)
            else:
                shap_abs = np.mean(np.abs(shap_values.reshape(shap_values.shape[0], -1)), axis=0)

            df = pd.DataFrame({'Gene': feature_names, 'Importance': shap_abs, 'Model': name})
            all_results.append(df)
            print(f"  ✓ Completed {name}")
        except Exception as ex:
            print(f"  ✗ Error with {name}: {ex}")

    if not all_results:
        raise RuntimeError("No SHAP results computed.")
    result = pd.concat(all_results, axis=0)
    return normalize_per_model(result)

def create_publication_plot(data, title, filename, plot_type="permutation"):
    mean_importance = data.groupby('Gene')['Importance'].mean().sort_values(ascending=False)
    top_genes = mean_importance.head(10).index.tolist()
    plot_data = data[data['Gene'].isin(top_genes)].copy()
    plot_data['Gene'] = pd.Categorical(plot_data['Gene'], categories=top_genes, ordered=True)
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    sns.barplot(data=plot_data, x='Gene', y='Importance', hue='Model',
                palette='Set2', ax=ax, edgecolor='black', linewidth=0.3)
    ax.set_title(title, fontsize=8, fontweight='bold', pad=10)
    ax.set_xlabel('Genes', fontsize=7, fontweight='bold')
    ax.set_ylabel(f'{plot_type.capitalize()} Importance', fontsize=7, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=6)
    ax.tick_params(axis='y', labelsize=6)
    legend = ax.legend(title='Model', fontsize=5, title_fontsize=6, frameon=True)
    legend.get_frame().set_linewidth(0.5)
    legend.get_frame().set_edgecolor('black')
    ax.grid(True, alpha=0.3, linewidth=0.3)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color('black')
    plt.tight_layout()
    plt.savefig(f"{filename}.svg", format='svg', bbox_inches='tight', pad_inches=0.05,
                facecolor='white', edgecolor='none')
    plt.show()
    return fig, ax

def create_summary_table(data, filename, analysis_type):
    pivot_table = data.pivot(index='Gene', columns='Model', values='Importance')
    pivot_table['Mean_Importance'] = pivot_table.mean(axis=1)
    pivot_table = pivot_table.sort_values('Mean_Importance', ascending=False)
    pivot_table = pivot_table.round(4)
    model_cols = [col for col in pivot_table.columns if col != 'Mean_Importance']
    pivot_table = pivot_table[model_cols + ['Mean_Importance']]
    pivot_table.to_csv(f"{filename}.csv")
    print(f"\n{analysis_type} Summary Table (Top 10):")
    print("="*50)
    print(pivot_table.head(10).to_string())
    return pivot_table

print("Starting Feature Importance Analysis (Permutation + SHAP, all models, normalized)...")
print("="*50)

try:
    print("\n1. Computing Permutation Importance all models...")
    perm_data = compute_permutation_all(best_models, X_train, y_train, feature_names,
                                        scoring='roc_auc', n_repeats=10, sample_size=1000)
    create_publication_plot(perm_data,
                            'Permutation Feature Importance (Top 10 Genes, Normalized)',
                            'perm_importance_allmodels', 'permutation')
    perm_table = create_summary_table(perm_data, 'perm_importance_allmodels', 'Permutation Importance')
except Exception as e:
    print(f"Permutation Importance Analysis failed: {e}")

try:
    print("\n2. Computing SHAP Importance all models...")
    shap_data = compute_shap_all(best_models, X_train, feature_names,
                                 sample_size=200, background_size=20)
    create_publication_plot(shap_data,
                            'SHAP Feature Importance (Top 10 Genes, Normalized)',
                            'shap_importance_allmodels', 'SHAP')
    shap_table = create_summary_table(shap_data, 'shap_importance_allmodels', 'SHAP Importance')
except Exception as e:
    print(f"SHAP Importance Analysis failed: {e}")

print("\nAnalysis completed!")
print("Files generated:")
print("- perm_importance_allmodels.svg")
print("- perm_importance_allmodels.csv")
print("- shap_importance_allmodels.svg")
print("- shap_importance_allmodels.csv")


###### END OF THE SCRIPT ######
