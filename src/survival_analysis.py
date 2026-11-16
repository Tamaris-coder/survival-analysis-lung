"""Survival analysis of the public NCCTG lung cancer dataset (lifelines)."""
import os
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.datasets import load_lung
from lifelines.statistics import logrank_test
from lifelines.utils import k_fold_cross_validation
from sklearn.model_selection import train_test_split

os.makedirs("results", exist_ok=True)

# 1. Load and prepare data
df = load_lung()
df.columns = df.columns.str.replace(".", "_", regex=False)
# lifelines may code status as 0/1 (1 = death) or 1/2 (2 = death)
df["event"] = (df["status"] == df["status"].max()).astype(int)
print("status values:", sorted(df["status"].unique()), "| sex values:", sorted(df["sex"].unique()))
df = df.drop(columns=["status", "inst"]).dropna()  # complete cases only
print("Patients:", len(df), "| Events:", df["event"].sum())

# 2. Kaplan-Meier, all patients
kmf = KaplanMeierFitter()
kmf.fit(df["time"], df["event"], label="All patients")
ax = kmf.plot_survival_function()
ax.set_xlabel("Days")
ax.set_ylabel("Survival probability")
plt.savefig("results/km_overall.png", dpi=200, bbox_inches="tight")
plt.close()

# 3. Kaplan-Meier by sex + log-rank test
fig, ax = plt.subplots()
for code, label in [(1, "Male"), (2, "Female")]:
    m = df["sex"] == code
    KaplanMeierFitter().fit(df.loc[m, "time"], df.loc[m, "event"],
                            label=label).plot_survival_function(ax=ax)
ax.set_xlabel("Days")
ax.set_ylabel("Survival probability")
plt.savefig("results/km_by_sex.png", dpi=200, bbox_inches="tight")
plt.close()

m = df["sex"] == 1
lr = logrank_test(df.loc[m, "time"], df.loc[~m, "time"],
                  df.loc[m, "event"], df.loc[~m, "event"])
print("Log-rank p-value (sex):", lr.p_value)

# 4. Cox proportional hazards model
cph = CoxPHFitter()
cph.fit(df, duration_col="time", event_col="event")
cph.print_summary()
cph.plot()
plt.savefig("results/cox_hazard_ratios.png", dpi=200, bbox_inches="tight")
plt.close()

# 5. Check the proportional hazards assumption
cph.check_assumptions(df, p_value_threshold=0.05, show_plots=False)

# 6. Train/test split: how well does the model rank patients? (C-index)
train, test = train_test_split(df, test_size=0.3, random_state=42)
cph_train = CoxPHFitter().fit(train, duration_col="time", event_col="event")
c_test = cph_train.score(test, scoring_method="concordance_index")
print("Test C-index:", round(c_test, 3))

# 7. Penalised (ridge) Cox with 5-fold cross-validation
for p in [0.0, 0.01, 0.1, 1.0]:
    scores = k_fold_cross_validation(
        CoxPHFitter(penalizer=p), df, duration_col="time",
        event_col="event", k=5, scoring_method="concordance_index", seed=42)
    print(f"penalizer={p}: mean CV C-index = {sum(scores)/len(scores):.3f}")