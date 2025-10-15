# Plotting Scripts

This directory contains visualization and plotting utilities for ICESEE results, with a focus on performance analysis and scaling studies.

## Scripts

### scaling_plots.py
Generate performance scaling plots for ICESEE parallel execution.

**Features**:
- Strong scaling analysis
- Weak scaling analysis
- Speedup calculations
- Efficiency plots
- Comparison across configurations

**Usage**:
```bash
python scaling_plots.py --data <scaling_data_file>
```

### scaling_plots_csv_details.py
Detailed scaling analysis using CSV data files.

**Features**:
- Read performance data from CSV
- Multiple metric visualization
- Detailed performance breakdowns
- Custom plot configurations
- Export publication-quality figures

**Usage**:
```bash
python scaling_plots_csv_details.py --csv <data.csv> --output <plots/>
```

## Data Format

Scaling scripts expect performance data in specific formats:
- Execution times
- Number of processors
- Problem sizes
- Memory usage
- Communication overhead

## Output

Plots are generated in common formats:
- PNG for quick viewing
- PDF for publications
- SVG for editing

## Use Cases

These scripts are particularly useful for:
- HPC performance analysis
- Parallel efficiency studies
- Configuration optimization
- Publication figures
- Technical reports

## Dependencies

- matplotlib (plotting)
- numpy (data processing)
- pandas (CSV handling)
- seaborn (enhanced visualizations, optional)
