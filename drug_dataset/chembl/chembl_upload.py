import os
import io
import json
import requests
import pandas as pd
from sqlalchemy import create_engine, text
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import psycopg2


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

MAX_WORKERS = 5        # threads for parallel API fetch
TABLE_NAME = "chembl_molecules"

# --- Connect to PostgreSQL ---
conn_str = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
engine = create_engine(conn_str)

# --- Create table schema if not exists ---
create_table_sql = f"""
DROP TABLE IF EXISTS {TABLE_NAME};

CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
    molecule_chembl_id TEXT PRIMARY KEY,
    molecule_data JSONB
);
"""

with engine.begin() as conn:
    conn.execute(text(create_table_sql))
    #conn.commit()


# --- COPY INSERT function ---
def copy_insert(engine, df, table_name):
    """Efficiently bulk insert a pandas DataFrame into PostgreSQL using COPY."""
    buffer = io.StringIO()
    df.to_csv(buffer, index=False, header=False, sep='\t')
    buffer.seek(0)

    conn = engine.raw_connection()
    cursor = conn.cursor()
    try:
        cursor.copy_expert(f"""
            COPY {table_name} (molecule_chembl_id, molecule_data)
            FROM STDIN WITH (FORMAT CSV, DELIMITER E'\t')
        """, buffer)
        conn.commit()
    except Exception as e:
        print("⚠️ COPY failed, rolling back:", e)
        conn.rollback()
        raise
    finally:
        cursor.close()
        conn.close()

# --- Fetch molecules in parallel ---
def fetch_and_insert(offset):
    """Fetch a batch from the ChEMBL API and insert into PostgreSQL."""
    params = {"limit": LIMIT, "offset": offset}
    resp = requests.get(API_URL, params=params)
    resp.raise_for_status()
    data = resp.json()

    if "molecules" not in data:
        print(f"Unexpected response at offset {offset}: {data}")
        return 0

    molecules = data["molecules"]
    df = pd.DataFrame([
        {"molecule_chembl_id": m["molecule_chembl_id"], "molecule_data": json.dumps(m)}
        for m in molecules
    ])
    copy_insert(engine, df, TABLE_NAME)
    return len(df)

# --- Main execution ---
print(f"🚀 Starting ChEMBL ingestion ({SAMPLE_SIZE:,} target molecules)...")
offsets = range(0, SAMPLE_SIZE, LIMIT)
inserted_total = 0

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = [executor.submit(fetch_and_insert, offset) for offset in offsets]

    with tqdm(total=SAMPLE_SIZE, desc="Downloading & inserting", unit="mol") as pbar:
        for f in as_completed(futures):
            try:
                inserted = f.result()
                inserted_total += inserted
                pbar.update(inserted)
            except Exception as e:
                print(f"⚠️ Error during batch insert: {e}")

print(f"✅ Done! {inserted_total:,} molecules successfully loaded into RDS.")