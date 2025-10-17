# train_permeability_model.py
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
import joblib
from tqdm import tqdm

from tdc.single_pred import ADME
data = ADME(name='Caco2_Wang')
df = data.get_data()  # returns a DataFrame with SMILES + measured permeability
# The dataset should have at least: "smiles" and "logPapp" (permeability)
print(df.head())
df = df.dropna()

df = df.rename(columns={
    "Drug": "smiles",
    "Y": "logPapp"
})

# 🧪 Compute RDKit descriptors
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "PSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotB": Descriptors.NumRotatableBonds(mol)
    }

desc_data = [compute_descriptors(smi) for smi in tqdm(df["smiles"], desc="Computing descriptors")]
desc_df = pd.DataFrame([d for d in desc_data if d is not None])
df = df.loc[desc_df.index]
df = pd.concat([df.reset_index(drop=True), desc_df.reset_index(drop=True)], axis=1)

# 📈 Train/test split
X = df[["MolWt", "LogP", "PSA", "HBD", "HBA", "RotB"]]
y = df["logPapp"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 🌲 Train model
model = RandomForestRegressor(n_estimators=300, random_state=42)
model.fit(X_train, y_train)

# 📊 Evaluate
preds = model.predict(X_test)
print(f"R² = {r2_score(y_test, preds):.3f}")
print(f"MAE = {mean_absolute_error(y_test, preds):.3f}")

# 💾 Save model
joblib.dump(model, "models/permeability_rf.joblib")
print("✅ Model saved to models/permeability_rf.joblib")
