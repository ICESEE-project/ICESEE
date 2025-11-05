"""
ICeSEE performance plotting utility.

Usage (inside a notebook or script):
    analyze_icesee("icesee_timing.log", scaling_type="strong_scaling")

Outputs:
  - PNGs and CSVs in icesee_perf_outputs
Notes:
  - Uses Assimilation time if present (preferred), else falls back to Wall-Clock,
    else to Computational Sum.
  - Baseline for speedup/efficiency is the smallest rank present (N0).
  - Strong scaling:
        speedup(N) = T(N0) / T(N)
        efficiency(N) = speedup(N) / (N/N0)
    Weak scaling (problem size ∝ N):
        efficiency(N) = T(N0) / T(N)
        speedup(N)    = (N/N0) * efficiency(N)
"""

import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ANSI = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')  # strip ANSI

def _time_to_seconds(s):
    """Convert 'DD:HH:MM:SS.mmm' to seconds (float), forgiving of odd chars."""
    s = s.strip()
    parts = s.split(":")
    if len(parts) == 4:
        d, h, m, rest = parts
    elif len(parts) == 3:        # tolerate 'HH:MM:SS.mmm'
        d, h, m, rest = 0, parts[0], parts[1], parts[2]
    else:
        d, h, m, rest = 0, 0, 0, parts[-1]
    def as_int(x):
        try: return int(x)
        except: return int(re.sub(r'[^0-9]','',str(x)) or 0)
    def as_float(x):
        try: return float(x)
        except:
            x = re.sub(r'[^0-9\.]', '', str(x))
            return float(x) if x else 0.0
    d,h,m,sec = as_int(d),as_int(h),as_int(m),as_float(rest)
    return (((d*24)+h)*60 + m)*60 + sec

def _time_to_minutes(s):
    return _time_to_seconds(s) / 60.0

def _time_to_hours(s):
    return _time_to_seconds(s) / 3600.0

def _find_time(window_lines, label_regex):
    for w in window_lines:
        if re.search(label_regex, w, re.I):
            tm = re.search(r"(\d{1,2}:\d{2}:\d{2}:\d{2}\.\d{3})", w)
            if not tm:
                tm = re.search(r"(\d+:\d{2}:\d{2}:\d{2}\.\d{3})", w)
            if tm:
                # return _time_to_seconds(tm.group(1))
                return _time_to_minutes(tm.group(1))
    return None

def _parse_log(log_path: str, model_np: int) -> pd.DataFrame:
    p = Path(log_path)
    text = p.read_text(errors="ignore")
    text = ANSI.sub("", text)  # remove ANSI escape codes
    lines = text.splitlines()

    # Find block headers that announce ranks
    start_idx = []
    for i, line in enumerate(lines):
        if "[ICESEE] Performance Metrics" in line or "[ICESEE] Metrics on" in line:
            start_idx.append(i)

    records = []
    for idx in start_idx:
        window = lines[idx: idx+60]
        header = window[0]
        ranks = None
        m = re.search(r"\((\d+)\s+ranks\)", header)
        if not m:
            m = re.search(r"Metrics on\s+(\d+)\s+ranks", header)
        if m:
            ranks = int(m.group(1))

        np = ranks/(model_np+1)
        ranks = int(np)

        wall = _find_time(window, r"Wall-Clock Time|Wall[- ]Clock Time")
        comp = _find_time(window, r"Computational Time")
        forecast = _find_time(window, r"Forecast Step Time")
        analysis = _find_time(window, r"Analysis Step Time")
        assim = _find_time(window, r"Assimilation Time")
        init_io = _find_time(window, r"Init file I/O Time")
        f_io   = _find_time(window, r"Forecast File I/O Time")
        a_io   = _find_time(window, r"Analysis File I/O Time")
        tot_io = _find_time(window, r"Total File I/O Time")
        f_noise= _find_time(window, r"Forecast Noise Time")
        f_mean = _find_time(window, r"Forecast Ensemble Mean Computation")
        a_mean = _find_time(window, r"Analysis Ensemble Mean Computation")

        if ranks is not None and (wall is not None or assim is not None or comp is not None):
            records.append({
                "ranks": ranks,
                "wallclock_s": wall,
                "computational_sum_s": comp,
                "assimilation_s": assim,
                "forecast_s": forecast or 0.0,
                "analysis_s": analysis or 0.0,
                "init_io_s": init_io or 0.0,
                "forecast_io_s": f_io or 0.0,
                "analysis_io_s": a_io or 0.0,
                "total_io_s": tot_io or 0.0,
                "forecast_noise_s": f_noise or 0.0,
                "forecast_mean_s": f_mean or 0.0,
                "analysis_mean_s": a_mean or 0.0,
            })

    df = pd.DataFrame(records).drop_duplicates()
    # If there are multiple entries per rank, keep the last one
    if not df.empty:
        df = df.sort_values("ranks").drop_duplicates(subset=["ranks"], keep="last").reset_index(drop=True)
    return df

def _compute_scaling(df: pd.DataFrame, scaling_type: str) -> pd.DataFrame:
    out = df.copy()

    # Prefer Assimilation time, else Wall-Clock, else Computational Sum
    metric_order = ["assimilation_s", "wallclock_s", "computational_sum_s"]
    for col in metric_order:
        if col in out.columns and out[col].notna().any():
            base_col = col
            break
    else:
        base_col = metric_order[-1]

    valid = out.dropna(subset=[base_col]).copy()
    if len(valid) < 2:
        out["speedup"] = float("nan")
        out["efficiency"] = float("nan")
        out["base_metric_col"] = base_col
        return out

    N0 = int(valid["ranks"].min())
    # N0=1
    T0 = float(valid.loc[valid["ranks"] == N0, base_col].iloc[0])

    out["baseline_ranks"] = N0
    out["baseline_time_s"] = T0
    out["metric_time_s"] = out[base_col]

    if scaling_type == "weak_scaling":
        out["efficiency"] = out["baseline_time_s"] / out["metric_time_s"]
        out["speedup"] = (out["ranks"] / N0) * out["efficiency"]
    else:  # strong scaling
        out["speedup"] = out["baseline_time_s"] / out["metric_time_s"]
        # out["efficiency"] = out["speedup"] / (out["ranks"] / N0)
        out["efficiency"] = (out["speedup"] / out["ranks"])*100
        # out["efficiency"] = (out["baseline_time_s"] / out["metric_time_s"]) / (out["ranks"] / N0) * 100

    out["base_metric_col"] = base_col
    return out

def _line_plot(x, y, xlabel, ylabel, title, out_path,flag=None):
    fig = plt.figure(figsize=(4.5, 4.5), dpi=300)
    # plt.plot(x, y, marker="o")
    if flag=='efficiency':
        ideal_efficiency = x/x
        plt.semilogx(x, y, 'o-b', markersize=6, linewidth=2, label='Measured')
        plt.semilogx(x, ideal_efficiency*100,  '--', markersize=8, linewidth=2, color='red',label='Ideal')
        plt.legend(['Measured','Ideal'], prop={'size': 10, 'weight': 'bold'})
    elif flag=='speedup':
        # npp = [2**i for i in range(0,len(x)+1)]
        # ideal_speedup = y[0]/npp
        # print(f"ideal_speedup: {ideal_speedup}, y: {y}, x: {x}")
        plt.loglog(x, y, 'o-b', markersize=6, linewidth=2, label='Measured')
        # plt.loglog(x, ideal_speedup, '--', markersize=8, linewidth=2, color='red', label='Ideal')
        # plt.legend(['Measured','Ideal'], prop={'size': 10, 'weight': 'bold'})
    else:
        plt.loglog(x, y, 'o-b', markersize=6, linewidth=2, label='Measured')
        npp = [2**i for i in range(0,len(x))]
        ideal_time= y[0]/npp
        # ideal_time= y[0]/x
        plt.loglog(x, ideal_time, '--', markersize=8, linewidth=2, color='red',label='Ideal')
        plt.legend(['Measured','Ideal'], prop={'size': 10, 'weight': 'bold'})
        
    plt.xlabel(xlabel, fontdict={'fontsize': 12, 'fontweight': 'bold'})
    plt.ylabel(ylabel, fontdict={'fontsize': 12, 'fontweight': 'bold'})
    plt.title(title, fontdict={'fontsize': 14, 'fontweight': 'bold'})
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.grid()
    # np = x.tolist()   
    pstr = ([f'{N:d}' for N in x])
    plt.xticks(x,pstr)
    if flag != 'efficiency' and flag != 'speedup':
        plt.yticks(y,[f'{N:.2f}' for N in y])

    if flag == 'speedup':
         ideal_speedup = y[0]/x
         plt.yticks(y,[f'{N:.2f}' for N in ideal_speedup])
    
    plt.gca().minorticks_off()

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def analyze_icesee(log_path: str, scaling_type: str = "strong_scaling", model_np: int = 1):
    out_dir = Path("icesee_perf_outputs")
    out_dir.mkdir(exist_ok=True, parents=True)

    df = _parse_log(log_path, model_np=model_np)  
    if df.empty:
        raise RuntimeError("No timing blocks found in the log. Check formatting or path.")

    metrics = _compute_scaling(df, scaling_type)

    # Save CSV
    csv_path = out_dir / f"summary_{scaling_type}.csv"
    metrics.to_csv(csv_path, index=False)

    # Plots: speedup & efficiency
    _line_plot(
        metrics["ranks"], metrics["speedup"],
        "Ranks (MPI processes)",
        f"Speedup (relative to N0={int(metrics['baseline_ranks'].iloc[0])})",
        f"{scaling_type.replace('_',' ').title()}: Speedup (metric: {metrics['base_metric_col'].iloc[0].replace('_s','')})",
        out_dir / f"{scaling_type}_speedup.png",
        flag='speedup'
    )
    _line_plot(
        metrics["ranks"], metrics["efficiency"],
        "Ranks (MPI processes)",
        "Parallel Efficiency",
        f"{scaling_type.replace('_',' ').title()}: Parallel Efficiency (metric: {metrics['base_metric_col'].iloc[0].replace('_s','')})",
        out_dir / f"{scaling_type}_efficiency.png",
        flag='efficiency'
    )

    # Per-phase timing (each on its own figure)
    if "forecast_s" in metrics:
        _line_plot(
            metrics["ranks"], metrics["forecast_s"],
            "Ranks (MPI processes)", "Time (mins)",
            "Forecast Step Time vs Ranks",
            out_dir / f"{scaling_type}_forecast_time.png"
        )
    if "analysis_s" in metrics:
        _line_plot(
            metrics["ranks"], metrics["analysis_s"],
            "Ranks (MPI processes)", "Time (mins)",
            "Analysis Step Time vs Ranks",
            out_dir / f"{scaling_type}_analysis_time.png"
        )
    # Prefer total I/O; if not present, plot forecast I/O
    io_series = None
    io_title = None
    if "total_io_s" in metrics and metrics["total_io_s"].notna().any():
        io_series = metrics["total_io_s"]
        io_title = "Total File I/O Time vs Ranks"
    elif "forecast_io_s" in metrics and metrics["forecast_io_s"].notna().any():
        io_series = metrics["forecast_io_s"]
        io_title = "Forecast File I/O Time vs Ranks"
    if io_series is not None:
        _line_plot(
            metrics["ranks"], io_series,
            "Ranks (MPI processes)", "Time (mins)",
            io_title,
            out_dir / f"{scaling_type}_io_time.png"
        )

    return {
        "csv": str(csv_path),
        "plots": [str(p) for p in sorted(out_dir.glob(f"{scaling_type}_*.png"))],
        "ranks": metrics["ranks"].tolist(),
        "baseline_ranks": int(metrics["baseline_ranks"].iloc[0]),
        "base_metric": metrics["base_metric_col"].iloc[0],
    }

# Example calls:
# analyze_icesee("icesee_timing.log", scaling_type="strong_scaling")
# analyze_icesee("icesee_timing.log", scaling_type="weak_scaling")