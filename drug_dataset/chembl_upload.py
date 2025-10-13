import os
print('hi')
import requests
import pandas as pd
from sqlalchemy import create_engine, text
from tqdm import tqdm
from jsonschema import validate
import json


'''
export DB_USER="postgres"
export DB_PASS="AzuleneLabs_2026"
export DB_HOST="azulene-1.cizeysmsgxmm.us-east-1.rds.amazonaws.com"
export DB_NAME="postgres"
export DB_PORT=5432
'''

# --- CONFIG ---
DB_USER = os.getenv("DB_USER")
DB_PASS = os.getenv("DB_PASS")
DB_HOST = os.getenv("DB_HOST")
DB_NAME = os.getenv("DB_NAME")
DB_PORT = os.getenv("DB_PORT", 5432)

API_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
BATCH_SIZE = 1000
TOTAL_RECORDS = 2000000        # total molecules in ChEMBL (approx)
SAMPLE_SIZE = int(TOTAL_RECORDS * 0.10)  # 10% sample (~200k)
LIMIT = BATCH_SIZE
OFFSET = 0

# --- Connect to PostgreSQL ---
conn_str = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
engine = create_engine(conn_str)

# --- Create table schema if not exists ---


create_table_sql = """
DROP TABLE IF EXISTS chembl_molecules;

CREATE TABLE IF NOT EXISTS chembl_molecules (
    molecule_chembl_id TEXT PRIMARY KEY,
    molecule_data JSONB

   
);
"""
with engine.begin() as conn:
    conn.execute(text(create_table_sql))
    #conn.commit()


# --- Fetch and insert data ---
print(f"Fetching {SAMPLE_SIZE} molecules from ChEMBL API...")

inserted = 0
pbar = tqdm(total=SAMPLE_SIZE, desc="Downloading & inserting")

while inserted < SAMPLE_SIZE:
    print("more hi")
    params = {"limit": LIMIT, "offset": OFFSET}
    resp = requests.get(API_URL, params=params)
    resp.raise_for_status()
    data = resp.json()

    if "molecules" not in data:
        print("Unexpected response")
        break
    molecules = data["molecules"]
    print(molecules[0]["molecule_chembl_id"])

    
    batch_size = 1000
    bi = 0

    with engine.begin() as conn:
        for i in range(0, len(molecules), batch_size):
            print(f"batch {bi}")
            bi += 1
            batch = molecules[i:i + batch_size]
            conn.execute(
                text("""
                    INSERT INTO chembl_molecules (molecule_chembl_id, molecule_data)
                    VALUES (:chembl_id, :data)
                    ON CONFLICT (molecule_chembl_id) DO NOTHING
                """),
                [
                    {"chembl_id": m["molecule_chembl_id"], "data": json.dumps(m)}
                    for m in batch
                ]
            )

    inserted += len(molecules)
    OFFSET += LIMIT
    pbar.update(len(molecules))

pbar.close()
print("✅ Done! Molecules loaded into RDS.")
