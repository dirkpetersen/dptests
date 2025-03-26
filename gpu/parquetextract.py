import json
import pandas as pd
from pathlib import Path

def process_json_data(row):
    """Extract and flatten nested GPU metrics from JSON data"""
    data = json.loads(row['json_data'])
    records = []
    
    for gpu in data.get('gpus', []):
        for process in gpu.get('processes', []):
            # Create base command without arguments
            base_cmd = process.get('command', '').split()[0] if process.get('command') else ''
            
            records.append({
                'hostname': data.get('hostname'),
                'timestamp': row['timestamp'],
                'gpuid': gpu.get('uuid'),
                'gputype': gpu.get('name'),
                'gpu_util': gpu.get('utilization.gpu'),
                'gpu_mem_used': gpu.get('memory.used'),
                'username': process.get('username'),
                'command': base_cmd,
                'pid': process.get('pid'),
                'gpu_memory_usage': process.get('gpu_memory_usage'),
                'job': f"{data['hostname']}-{process.get('username')}-{base_cmd}-{process.get('pid')}"
            })
    
    return records

def main(input_parquet, output_csv):
    df = pd.read_parquet(input_parquet)
    all_records = []
    
    for _, row in df.iterrows():
        all_records.extend(process_json_data(row))
    
    output_df = pd.DataFrame(all_records)
    output_df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Input parquet file')
    parser.add_argument('output', help='Output CSV file')
    args = parser.parse_args()
    
    main(args.input, args.output)
