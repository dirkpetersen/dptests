import pandas as pd
from datetime import datetime

def summarize_gpu_stats(input_csv, output_csv):
    # Read CSV with proper datetime parsing
    df = pd.read_csv(input_csv, parse_dates=['timestamp'])
    
    # Add column verification and normalization
    original_columns_set = set(df.columns)
    column_normalized = False
    
    if 'Job' not in df.columns:  
        df.columns = df.columns.str.strip().str.lower()
        column_normalized = True
        if 'job' not in df.columns:
            raise ValueError("CSV missing required 'Job/job column")

    # Determine actual column names being used
    job_col = 'job'  # After normalization or if original was lowercase
    timestamp_col = 'timestamp'
    
    # Define aggregation functions using verified column names
    agg_dict = {
        timestamp_col: ['min', 'max'],
        'gpu_util': 'mean',
        'gpu_mem_used': 'mean',
        'gpu_memory_usage': 'mean',
        'hostname': 'first',     # Changed from 'machine'
        'username': 'first',     # Changed from 'user'
        'command': 'first',      # Add command
        **({'gputype': 'first'} if 'gputype' in df.columns else {})  # Handle alternative names
    }

    # Group by the correct job column name
    grouped = df.groupby(job_col).agg(agg_dict)
    
    # Flatten multi-index columns and rename time columns
    grouped.columns = ['_'.join(col).strip('_') for col in grouped.columns.values]
    
    # Rename timestamp columns FIRST
    grouped = grouped.rename(columns={
        f'{timestamp_col}_min': 'start_time',
        f'{timestamp_col}_max': 'end_time',
        'hostname_first': 'hostname',    # Clean up column names
        'username_first': 'username',
        'command_first': 'command'
    })
    
    # NOW define column order with FINAL names
    column_order = [
        'start_time', 'end_time',
        'duration_min',
        'gpu_util_mean', 'gpu_mem_used_mean', 'gpu_memory_usage_mean',
        'hostname', 'username', 'command'
    ]
    if 'gputype_first' in grouped.columns:
        column_order.append('gputype_first')

    # Calculate duration in minutes
    grouped['duration_min'] = (grouped['end_time'] - grouped['start_time']).dt.total_seconds() / 60
    
    grouped = grouped[column_order].reset_index()
    
    # Save to CSV
    grouped.to_csv(output_csv, index=False, float_format='%.2f')

if __name__ == '__main__':
    summarize_gpu_stats('gpu_stats.csv', 'gpu_job_summary.csv')
