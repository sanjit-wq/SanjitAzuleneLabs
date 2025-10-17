import os
import json
import requests
import pandas as pd
from sqlalchemy import create_engine, text
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import joblib

# --- CONFIG ---
DB_USER = os.getenv("DB_USER")
DB_PASS = os.getenv("DB_PASS")
DB_HOST = os.getenv("DB_HOST")
DB_NAME = os.getenv("DB_NAME")
DB_PORT = os.getenv("DB_PORT", 5432)

BATCH_SIZE = 1000

permeability_model = joblib.load("models/permeability_rf.joblib")


# --- Connect to PostgreSQL ---
conn_str = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
engine = create_engine(conn_str)


# --- RDKit-based calculations ---
def compute_rdkit_features(smiles):
    """Compute molecular weight, logP, etc. using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        return {
            "mol_weight": Descriptors.MolWt(mol),
            "logp_rdkit": Crippen.MolLogP(mol)
        }
    except Exception:
        return None


# --- PubChem fallback API ---
def fetch_pubchem_logp(smiles):
    """Fetch experimental logP from PubChem."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/LogP,IsomericSMILES/JSON"
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            props = resp.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            return props.get("LogP")
    except Exception:
        return None


# --- Load data ---
def fetch_data(limit=None):
    """Fetch data from PostgreSQL table."""
    query = "SELECT * FROM drug_properties"
    if limit:
        query += f" LIMIT {limit}"
    return pd.read_sql(query, engine)


# --- Enrichment Pipeline ---
def enrich_dataframe(df):

    df = df.dropna(subset=["smiles"]).copy()
    df = df[df["smiles"].apply(lambda s: isinstance(s, str) and len(s.strip()) > 0)]
    
    enriched = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Enriching molecules"):
        smiles = row.get("smiles")

        mol = Chem.MolFromSmiles(smiles)
        logp = Descriptors.MolLogP(mol)
        psa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        desc = [
            mw,
            logp,
            psa,
            hbd,
            hba,
            Descriptors.NumRotatableBonds(mol)
    ]
        permeability = float(permeability_model.predict([desc])[0])

        row["permeability"] = permeability
        row["hba"]  = hba
        row["hbd"]  = hbd
        row["psa"]  = psa
        row["molecular_weight"] = mw



        data_origin = {}

        if not smiles:
            enriched.append(row)
            continue

        # --- Compute with RDKit ---
        rdkit_features = compute_rdkit_features(smiles)
        if rdkit_features:
            if pd.isna(row.get("logp")) and rdkit_features["logp_rdkit"] is not None:
                row["logp"] = rdkit_features["logp_rdkit"]
                data_origin["logp"] = "rdkit"
            if pd.isna(row.get("binding_free_energy")):
                # (placeholder example)
                row["binding_free_energy"] = -0.1 * rdkit_features["logp_rdkit"]
                data_origin["binding_free_energy"] = "estimated_rdkit"

        # --- PubChem fallback ---
        if pd.isna(row.get("logp")):
            logp_pubchem = fetch_pubchem_logp(smiles)
            if logp_pubchem is not None:
                row["logp"] = logp_pubchem
                data_origin["logp"] = "pubchem"

        # --- Simple pKa estimation (toy model) ---
        if pd.isna(row.get("pka")) and rdkit_features:
            row["pka"] = 7.0 - 0.2 * rdkit_features["logp_rdkit"]
            data_origin["pka"] = "estimated_rdkit"

        # --- Solubility fallback (basic logS estimation) ---
        if pd.isna(row.get("solubility")) and rdkit_features:
            logp = rdkit_features["logp_rdkit"]
            molwt = rdkit_features["mol_weight"]
            row["solubility"] = -0.01 * molwt - 0.5 * logp
            data_origin["solubility"] = "estimated_rdkit"

        # Track origin of each field
        row["metadata"] = json.dumps({"data_origin": data_origin})
        enriched.append(row)

    return pd.DataFrame(enriched)


# --- Save enriched data ---
def save_to_postgres(df):
    ## Make sure all metadata rows are strings
    df["metadata"] = df["metadata"].apply(
        lambda x: json.dumps(x) if isinstance(x, dict) else x
    )
    
    df.to_sql("drug_properties_enriched", engine, if_exists="replace", index=False)
    print(f"✅ Saved {len(df)} enriched rows to drug_properties_enriched table.")


# --- Main workflow ---
if __name__ == "__main__":
    print("🚀 Loading data from database...")
    df = fetch_data(limit=5000)  # adjust for testing
    print(f"Loaded {len(df)} molecules.")

    print("🔬 Enriching data...")
    enriched_df = enrich_dataframe(df)

    print("💾 Saving to database...")
    save_to_postgres(enriched_df)

    print("✅ Enrichment complete!")
