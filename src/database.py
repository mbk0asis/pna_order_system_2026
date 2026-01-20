import sqlite3
import pandas as pd
from datetime import datetime
import os

DB_PATH = "pna_design.db"

def get_connection():
    """Returns a connection to the SQLite database."""
    return sqlite3.connect(DB_PATH)

def init_db():
    """Initialize the database tables."""
    conn = get_connection()
    c = conn.cursor()
    
    # Genes Table
    c.execute('''
        CREATE TABLE IF NOT EXISTS genes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            symbol TEXT UNIQUE NOT NULL,
            accession TEXT,
            description TEXT,
            sequence TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Probes Table (Results)
    c.execute('''
        CREATE TABLE IF NOT EXISTS probes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT,
            sequence TEXT,
            start_pos INTEGER,
            length INTEGER,
            tm_giesen REAL,
            gc_content REAL,
            purine_content REAL,
            valid BOOLEAN,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(gene_symbol) REFERENCES genes(symbol)
        )
    ''')
    
    conn.commit()
    conn.close()

def save_gene(symbol, accession, description, sequence):
    """Saves or updates a gene record."""
    conn = get_connection()
    c = conn.cursor()
    try:
        c.execute('''
            INSERT INTO genes (symbol, accession, description, sequence)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(symbol) DO UPDATE SET
            accession=excluded.accession,
            description=excluded.description,
            sequence=excluded.sequence,
            created_at=CURRENT_TIMESTAMP
        ''', (symbol.upper(), accession, description, sequence))
        conn.commit()
    except Exception as e:
        print(f"DB Error save_gene: {e}")
    finally:
        conn.close()

def get_gene(symbol):
    """Retrieves gene details if it exists."""
    conn = get_connection()
    c = conn.cursor()
    c.execute("SELECT sequence, accession, description FROM genes WHERE symbol = ?", (symbol.upper(),))
    row = c.fetchone()
    conn.close()
    return row if row else None

def save_probes(gene_symbol, df_probes):
    """Saves a dataframe of designed probes to the database."""
    if df_probes.empty:
        return
        
    conn = get_connection()
    
    # Prepare data for insertion
    # We want to keep a history, so we just insert. 
    # Optionally we could clear old probes for this gene first if we want only latest?
    # For now, let's keep history but maybe add a batch ID? 
    # Simpler: just insert.
    
    # Select relevant columns that match our schema
    # Rename columns to match schema if necessary, or just map them
    data = df_probes.copy()
    data['gene_symbol'] = gene_symbol.upper()
    data['created_at'] = datetime.now()
    
    # Map DataFrame columns to DB columns
    # DB: sequence, start_pos, length, tm_giesen, gc_content, purine_content, valid
    # DF: probe_sequence, start, length, tm_pna_giesen, gc_content, purine_content, valid
    
    db_data = data[['gene_symbol', 'probe_sequence', 'start', 'length', 'tm_pna_giesen', 
                    'gc_content', 'purine_content', 'valid', 'created_at']]
    
    db_data.columns = ['gene_symbol', 'sequence', 'start_pos', 'length', 'tm_giesen', 
                       'gc_content', 'purine_content', 'valid', 'created_at']
    
    db_data.to_sql('probes', conn, if_exists='append', index=False)
    conn.close()

def get_all_genes():
    """Returns a list of all saved genes."""
    conn = get_connection()
    c = conn.cursor()
    c.execute("SELECT symbol, accession, created_at FROM genes ORDER BY created_at DESC")
    rows = c.fetchall()
    conn.close()
    return rows
