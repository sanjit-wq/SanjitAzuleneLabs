"""
merge_datasets.py (updated for ChEMBL nested molecular properties)
"""

import os
import json
import requests
import pandas as pd
from sqlalchemy import create_engine, text
from tqdm import tqdm


# ----------------------------
# Environment Config
# ----------------------------
DB_USER = os.getenv("DB_USER")
DB_PASS = os.getenv("DB_PASS")
DB_HOST = os.getenv("DB_HOST")
DB_NAME = os.getenv("DB_NAME")
DB_PORT = os.getenv("DB_PORT", 5432)

# ----------------------------
# Database Connection
# ----------------------------
conn_str = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
engine = create_engine(conn_str)

# ----------------------------
# Schema Setup
# ----------------------------
create_table_sql = """
DROP TABLE IF EXISTS drug_properties;

CREATE TABLE IF NOT EXISTS drug_properties (
    chembl_id TEXT,
    smiles TEXT,
    binding_free_energy FLOAT,
    solubility FLOAT,
    logp FLOAT,
    permeability FLOAT,
    pka FLOAT,
    molecular_weight FLOAT,
    hba INT,
    hbd INT,
    psa FLOAT,
    rtb INT,
    qed_weighted FLOAT,
    source TEXT,
    metadata JSONB,
    PRIMARY KEY (chembl_id, source)
);
"""
with engine.begin() as conn:
    conn.execute(text(create_table_sql))
print("✅ Database schema initialized.")


# ----------------------------
# Safe Insert
# ----------------------------
def insert_into_db(df, source_name, batch_size=500):
    df = df.copy()
    df["source"] = source_name
    df = df.dropna(subset=["chembl_id"])

    if "metadata" in df.columns:
        df["metadata"] = df["metadata"].apply(
            lambda x: json.dumps(x) if isinstance(x, (dict, list)) else x
        )

    for start in range(0, len(df), batch_size):
        end = start + batch_size
        chunk = df.iloc[start:end]
        try:
            chunk.to_sql("drug_properties", engine, if_exists="append", index=False, method="multi")
        except Exception as e:
            print(f"⚠️ Batch insert failed ({start}-{end}): {e}")


# ----------------------------
# Fetch ChEMBL Data (updated)
# ----------------------------
def fetch_chembl_properties(limit=1000, total=20000):
    api_url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    dfs = []
    offset = 0
    pbar = tqdm(total=total, desc="Fetching ChEMBL molecules")

    while offset < total:
        resp = requests.get(api_url, params={"limit": limit, "offset": offset})
        resp.raise_for_status()
        data = resp.json().get("molecules", [])
        if not data:
            break

        df = pd.json_normalize(data)

        # Flatten nested molecule_properties
        if "molecule_properties" in df.columns:
            props = pd.json_normalize(df["molecule_properties"])
            props.columns = [c.split(".")[-1] for c in props.columns]
            df = pd.concat([df, props], axis=1)

        # Select relevant fields
        selected_cols = {
            "molecule_chembl_id": "chembl_id",
            "molecule_structures.canonical_smiles": "smiles",
            "alogp": "logp",
            "full_mwt": "molecular_weight",
            "cx_most_apka": "pka",
            "psa": "psa",
            "hba": "hba",
            "hbd": "hbd",
            "rtb": "rtb",
            "qed_weighted": "qed_weighted",
        }

        # Keep available columns only
        selected_cols = {k: v for k, v in selected_cols.items() if k in df.columns}

        subset = df[list(selected_cols.keys())].rename(columns=selected_cols)

        # Extract some metadata (formula, MW, etc.)
        subset["metadata"] = df.apply(
            lambda x: {
                "full_molformula": x.get("full_molformula"),
                "aromatic_rings": x.get("aromatic_rings"),
                "num_ro5_violations": x.get("num_ro5_violations"),
                "heavy_atoms": x.get("heavy_atoms"),
            },
            axis=1
        )

        # Add missing numeric placeholders
        for col in ["binding_free_energy", "solubility", "permeability"]:
            subset[col] = None

        dfs.append(subset)
        offset += limit
        pbar.update(limit)

    pbar.close()
    chembl_df = pd.concat(dfs, ignore_index=True)
    insert_into_db(chembl_df, "ChEMBL")
    print(f"✅ Inserted {len(chembl_df)} ChEMBL molecules.")


# ----------------------------
# Fetch AqSolDB Data
# ----------------------------
def fetch_aqsol_properties(file_path="AqSol/AqSolDB.csv"):
    """
    Load AqSolDB, clean and standardize columns, and insert into SQL.
    Designed to align with the unified `drug_properties` schema.
    """

    print("📦 Loading AqSolDB dataset...")
    aqsol = pd.read_csv(file_path)

    # --- Rename key columns for consistency ---
    aqsol = aqsol.rename(columns={
        "SMILES": "smiles",
        "Solubility": "solubility",
        "MolWt": "molecular_weight",
        "MolLogP": "logp",
        "NumHAcceptors": "hba",
        "NumHDonors": "hbd",
        "TPSA": "psa",
        "RingCount": "rtb"  # approximate mapping; true 'rotatable bonds' exist as NumRotatableBonds
    })

    # --- Use more direct rotatable bond field if available ---
    aqsol["rtb"] = aqsol["NumRotatableBonds"]

    # --- Generate synthetic chembl-like ID ---
    aqsol["chembl_id"] = aqsol["InChIKey"].fillna("").apply(
        lambda x: f"AQSOL_{x}" if x else f"AQSOL_{abs(hash(x)) % 10**10}"
    )

    # --- Create optional placeholders for other target properties ---
    optional_cols = ["binding_free_energy", "permeability", "pka"]
    for col in optional_cols:
        aqsol[col] = None

    # --- Construct metadata JSON ---
    meta_cols = [
        "Name", "InChI", "Occurrences", "Group", "HeavyAtomCount", "NumHeteroatoms",
        "NumValenceElectrons", "NumAromaticRings", "NumSaturatedRings",
        "NumAliphaticRings", "LabuteASA", "BalabanJ", "BertzCT"
    ]

    def make_metadata(row):
        meta = {col: row.get(col) for col in meta_cols if col in row}
        meta["source"] = "AqSolDB"
        return json.dumps(meta)

    aqsol["metadata"] = aqsol.apply(make_metadata, axis=1)

    # --- Select and reorder final columns for database consistency ---
    final_cols = [
        "chembl_id", "smiles", "binding_free_energy", "solubility", "logp",
        "permeability", "pka", "molecular_weight", "hba", "hbd", "psa", "rtb", "metadata"
    ]

    aqsol_final = aqsol[final_cols].drop_duplicates(subset=["chembl_id", "smiles"])

    print(f"🚀 Inserting {len(aqsol_final)} AqSolDB entries into database...")
    insert_into_db(aqsol_final, "AqSolDB")
    print("✅ AqSolDB data successfully merged.")



# ----------------------------
# Main Routine
# ----------------------------
if __name__ == "__main__":
    print("🚀 Starting data merge pipeline...")
    fetch_chembl_properties(limit=1000, total=20000)
    fetch_aqsol_properties("AqSol/AqSolDB.csv")
    print("🎉 All datasets merged into 'drug_properties' table.")
