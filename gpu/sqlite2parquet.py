import argparse
import os
import duckdb

def find_sqlite_files(folder_path):
    sqlite_files = []
    for root, _, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.sqlite'):
                sqlite_files.append(os.path.join(root, file))
    return sqlite_files

def has_gpu_stats(conn, sqlite_path):
    query = "SELECT COUNT(*) FROM SQLITE_SCAN(?, 'sqlite_master') WHERE type='table' AND name='gpu_stats'"
    return conn.execute(query, [sqlite_path]).fetchone()[0] > 0

def process_files(sqlite_paths, conn):
    merged_created = False
    for path in sqlite_paths:
        if not has_gpu_stats(conn, path):
            continue
        if not merged_created:
            conn.execute("CREATE TABLE merged_stats AS SELECT * FROM SQLITE_SCAN(?, 'gpu_stats')", [path])
            merged_created = True
        else:
            conn.execute("INSERT INTO merged_stats SELECT * FROM SQLITE_SCAN(?, 'gpu_stats')", [path])
    return merged_created

def main():
    parser = argparse.ArgumentParser(description='Merge GPU stats from SQLite files to Parquet')
    parser.add_argument('folder', help='Path to folder containing SQLite files')
    args = parser.parse_args()
    
    if not os.path.isdir(args.folder):
        raise ValueError(f"Invalid folder path: {args.folder}")
    
    with duckdb.connect() as conn:
        sqlite_files = find_sqlite_files(args.folder)
        if not sqlite_files:
            print("No SQLite files found")
            return

        has_data = process_files(sqlite_files, conn)
        
        if has_data:
            conn.execute("COPY merged_stats TO 'gpu_stats_merged.parquet' (FORMAT PARQUET)")
            print("Merged Parquet file created: gpu_stats_merged.parquet")
        else:
            print("No valid gpu_stats tables found")

if __name__ == '__main__':
    main()
