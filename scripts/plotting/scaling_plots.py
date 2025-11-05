import re
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict

def parse_time(time_str: str) -> float:
    """Parse time strings in DAY:HR:MIN:SEC.ms, MIN:SEC.ms, or SEC.ms to seconds."""
    if '.' not in time_str:
        return 0.0
    parts = time_str.split('.')
    ms = int(parts[1][:3]) / 1000  # Take up to 3 digits for ms
    time_parts = parts[0].split(':')
    
    # Initialize components
    d, h, m, s = 0, 0, 0, 0
    
    # Assign values based on number of components
    if len(time_parts) == 4:
        d, h, m, s = map(int, time_parts)
    elif len(time_parts) == 3:
        h, m, s = map(int, time_parts)
    elif len(time_parts) == 2:
        m, s = map(int, time_parts)
    elif len(time_parts) == 1:
        s = int(time_parts[0])
    else:
        return 0.0
    
    total_sec = d * 86400 + h * 3600 + m * 60 + s + ms
    return total_sec

def extract_metrics_from_log(filename: str) -> pd.DataFrame:
    """Parse the log file and extract metrics for blocks with rank information."""
    with open(filename, 'r') as f:
        content = f.read()
    
    # Remove ANSI escape codes
    content = re.sub(r'\x1b\[[0-9;]*m', '', content)
    
    # Split into lines
    lines = content.splitlines()
    
    data: List[Dict] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('╔═') and i + 1 < len(lines) and '[ICESEE] Performance Metrics' in lines[i+1]:
            # Start of a block
            block_lines = []
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith('╚═'):
                block_lines.append(lines[j].strip())
                j += 1
            i = j
            
            # Process block
            title_line = block_lines[0][2:-2].strip() if block_lines else ''  # Remove ║ and spaces
            ranks = None
            match = re.search(r'\((\d+) ranks\)', title_line)
            if match:
                ranks = int(match.group(1))
            
            if ranks is None:
                continue  # Skip blocks without rank info
            
            metrics = {}
            for bline in block_lines[1:]:
                if bline.startswith('║') and ':' in bline:
                    parts = bline[2:].split(':', 1)  # Remove ║ 
                    if len(parts) == 2:
                        metric_name = parts[0].strip()
                        value_part = parts[1].strip().split(' ')[0]  # Get the time string
                        metrics[metric_name] = parse_time(value_part)
            
            if metrics:
                entry = {
                    'ranks': ranks,
                    'wall_sec': metrics.get('Wall-Clock Time (max)', 0.0),
                    'forecast_sec': metrics.get('Forecast Step Time', 0.0),
                    'analysis_sec': metrics.get('Analysis Step Time', 0.0),
                    'io_sec': metrics.get('Total File I/O Time', 0.0)
                }
                data.append(entry)
        else:
            i += 1
    
    df = pd.DataFrame(data)
    df = df.sort_values(by='ranks').reset_index(drop=True)
    return df

def plot_metrics(df: pd.DataFrame, scaling_type: str):
    """Plot the metrics based on scaling type."""
    if df.empty:
        print("No data to plot.")
        return
    
    ref_p = df['ranks'].min()
    ref_t = df.loc[df['ranks'] == ref_p, 'wall_sec'].values[0]
    
    if scaling_type == 'strong_scaling':
        df['speedup'] = ref_t / df['wall_sec']
        df['efficiency'] = (ref_t / df['wall_sec']) / (df['ranks'] / ref_p) * 100
    elif scaling_type == 'weak_scaling':
        df['speedup'] = df['ranks'] / ref_p * (ref_t / df['wall_sec'])
        df['efficiency'] = ref_t / df['wall_sec'] * 100
    else:
        raise ValueError("Invalid scaling_type. Choose 'strong_scaling' or 'weak_scaling'.")
    
    # Plot 1: Wall-Clock Time vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['wall_sec'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Wall-Clock Time (seconds)')
    plt.title(f'Wall-Clock Time vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()
    
    # Plot 2: Speedup vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['speedup'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Speedup')
    plt.title(f'Speedup vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()
    
    # Plot 3: Parallel Efficiency vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['efficiency'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Parallel Efficiency (%)')
    plt.title(f'Parallel Efficiency vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()
    
    # Plot 4: Forecast Step Time vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['forecast_sec'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Forecast Step Time (seconds)')
    plt.title(f'Forecast Step Time vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()
    
    # Plot 5: Analysis Step Time vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['analysis_sec'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Analysis Step Time (seconds)')
    plt.title(f'Analysis Step Time vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()
    
    # Plot 6: Total File I/O Time vs Ranks
    plt.figure(figsize=(10, 6))
    plt.plot(df['ranks'], df['io_sec'], marker='o')
    plt.xlabel('Number of Ranks')
    plt.ylabel('Total File I/O Time (seconds)')
    plt.title(f'Total File I/O Time vs Number of Ranks ({scaling_type})')
    plt.grid(True)
    plt.show()

# Example usage:
# scaling_type = 'weak_scaling'  # or 'strong_scaling'
# df = extract_metrics_from_log('icesee_timing.log')
# plot_metrics(df, scaling_type)