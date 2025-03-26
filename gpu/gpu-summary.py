import pandas as pd
from datetime import datetime

def summarize_gpu_stats(input_csv, output_csv):
    # Read CSV with proper datetime parsing
    df = pd.read_csv(input_csv, parse_dates=['timestamp'])
    
    # Define aggregation functions for columns
    agg_dict = {
        'timestamp': ['min', 'max'],  # Start and end times
        'gpu_util': 'mean',
        'gpu_mem_used': 'mean',
        'gpu_memory_usage': 'mean',
        # Preserve other columns (take first occurrence assuming they're consistent)
        'machine': 'first',
        'user': 'first',
        'gpu_model': 'first'
    }
    
    # Group by Job and aggregate
    grouped = df.groupby('Job').agg(agg_dict)
    
    # Flatten multi-index columns and rename time columns
    grouped.columns = ['_'.join(col).strip('_') for col in grouped.columns.values]
    grouped = grouped.rename(columns={
        'timestamp_min': 'start_time',
        'timestamp_max': 'end_time'
    })
    
    # Calculate duration in minutes
    grouped['duration_min'] = (grouped['end_time'] - grouped['start_time']).dt.total_seconds() / 60
    
    # Reorder columns
    column_order = [
        'start_time', 'end_time', 'duration_min',
        'gpu_util_mean', 'gpu_mem_used_mean', 'gpu_memory_usage_mean',
        'machine_first', 'user_first', 'gpu_model_first'
    ]
    grouped = grouped[column_order].reset_index()
    
    # Save to CSV
    grouped.to_csv(output_csv, index=False, float_format='%.2f')

if __name__ == '__main__':
    summarize_gpu_stats('gpu_stats.csv', 'gpu_job_summary.csv')
