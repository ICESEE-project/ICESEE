# author: brian kyanjo
# date: september 2025
# This script processes ICESEE log files to extract performance metrics and generate scaling plots.
# - Only compute speedup/efficiency for WALLCLOCK and FORECAST metrics
# - Add stacked bar plots for (Analysis, Forecast, Ensemble Init, Assimilation, Total I/O)
#
# It writes PNGs/CSVs under icesee_perf_outputs and returns a small manifest.

import re
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ANSI = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')  # strip ANSI

def _time_to_seconds(s):
    s = s.strip()
    parts = s.split(":")
    if len(parts) == 4:
        d, h, m, rest = parts
    elif len(parts) == 3:
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

def _find_time(window_lines, label_regex):
    for w in window_lines:
        if re.search(label_regex, w, re.I):
            tm = re.search(r"(\d{1,2}:\d{2}:\d{2}:\d{2}\.\d{3})", w)
            if not tm:
                tm = re.search(r"(\d+:\d{2}:\d{2}:\d{2}\.\d{3})", w)
            if tm:
                return _time_to_minutes(tm.group(1))
    return None

def parse_icesee_log(log_path: str, model_np: int = 1) -> pd.DataFrame:
    p = Path(log_path)
    text = p.read_text(errors="ignore")
    text = ANSI.sub("", text)
    lines = text.splitlines()

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

        # adjust ranks to number of model subprocs if user used nested mpi (optional)
        np_eff = ranks/(model_np+1)
        ranks_eff = int(np_eff)

        wall = _find_time(window, r"Wall-Clock Time|Wall[- ]Clock Time")
        comp = _find_time(window, r"Computational Time")
        forecast = _find_time(window, r"Forecast Step Time")
        analysis = _find_time(window, r"Analysis Step Time")
        assim = _find_time(window, r"Assimilation Time")
        init_io = _find_time(window, r"Ensemble Init Time")
        tot_io = _find_time(window, r"Total File I/O Time")
        init_io = _find_time(window, r"Init file I/O Time")
        f_io   = _find_time(window, r"Forecast File I/O Time")
        a_io   = _find_time(window, r"Analysis File I/O Time")

        if ranks is not None and (wall is not None or assim is not None or comp is not None):
            records.append({
                "ranks": ranks_eff,
                "wall_m": wall,
                "forecast_m": forecast or 0.0,
                "analysis_m": analysis or 0.0,
                "assim_m": assim,
                "ens_init_m": init_io or 0.0,
                "io_m": tot_io or 0.0,
                "comp_sum_m": comp,
                "io_i_m": init_io,
                "io_f_m": f_io,
                "io_a_m": a_io
            })

    df = pd.DataFrame(records).drop_duplicates()
    if not df.empty:
        df = df.sort_values("ranks").drop_duplicates(subset=["ranks"], keep="last").reset_index(drop=True)
    return df

def compute_scaling(df: pd.DataFrame, scaling_type: str, column: str) -> pd.DataFrame:
    """Compute speedup/efficiency for a specific timing column (minutes)."""
    out = df[["ranks", column]].dropna().copy()
    out.rename(columns={column: "time_m"}, inplace=True)
    if out.empty or len(out) < 2:
        out["speedup"] = float("nan")
        out["efficiency"] = float("nan")
        return out

    N0 = int(out["ranks"].min())
    T0 = float(out.loc[out["ranks"] == N0, "time_m"].iloc[0])
    out["baseline_ranks"] = N0
    out["baseline_time_m"] = T0

    if scaling_type == "weak_scaling":
        out["efficiency"] = (T0 / out["time_m"]) * 100.0
        out["speedup"] = (out["ranks"] / N0) * (out["efficiency"] / 100.0)
    else:
        speedup = T0 / out["time_m"]
        out["speedup"] = speedup
        out["efficiency"] = (speedup / (out["ranks"] / N0)) * 100.0
    return out

def plot_line(x, y, xlabel, ylabel, title, out_path, ideal=None, logx=True, logy=False, _loglog=False, model_np=1):
    fig = plt.figure(figsize=(4.5, 4.5), dpi=300)
    np=x
    x = x*(model_np+1)  # convert back to total ranks for plotting
    if _loglog:
        plt.loglog(x, y, marker="o", markersize=6, linewidth=2, color='b',label="Measured")
        if ideal is not None:
            # ideal = (x / x[0]) * y[0]
            plt.loglog(x, ideal, linestyle="--", markersize=6, linewidth=2, color='r',label="Ideal")
            plt.legend(['Measured','Ideal'], prop={'size': 10, 'weight': 'bold'})
            # speed_up = np.arange(1, len(x)+1)
            # plt.yticks(y,[f'{N:.2f}' for N in y])
            # over ride yticks to show speedup values
            # plt.yticks(y,[f'{N:d}' for N in y])
    elif logx:
        plt.semilogx(x, y, marker="o", markersize=6, linewidth=2, color='b',label="Measured")
        if ideal is not None:
            plt.semilogx(x, ideal, linestyle="--",markersize=6, linewidth=2, color='r',label="Ideal")
            plt.legend(['Measured','Ideal'], prop={'size': 10, 'weight': 'bold'})
    else:
        if logy:
            plt.semilogy(x, y, marker="o", markersize=6, linewidth=2, color='b',label="Measured")
            if ideal is not None:
                plt.semilogy(x, ideal, linestyle="--", markersize=6, linewidth=2, color='r',label="Ideal")
                plt.legend(prop={'size': 10, 'weight': 'bold'})
        else:
            plt.plot(x, y, marker="o", markersize=6, linewidth=2, color='b',label="Measured")
            if ideal is not None:
                plt.plot(x, ideal, linestyle="--", markersize=6, linewidth=2, color='r',label="Ideal")
                plt.legend(prop={'size': 10, 'weight': 'bold'})
    plt.xlabel(xlabel, fontdict={'fontsize': 12, 'fontweight': 'bold'})
    plt.ylabel(ylabel, fontdict={'fontsize': 12, 'fontweight': 'bold'})
    plt.title(title, fontdict={'fontsize': 12, 'fontweight': 'bold'})
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.grid(True, which="both", axis="both", linestyle="--", alpha=0.5)
    pstr = ([f'{N:d}' for N in x])
    plt.xticks(x,pstr)
    plt.gca().minorticks_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=200)
    plt.close(fig)

def plot_stacked_bars_(df: pd.DataFrame, title: str, out_path: Path, model_np=1):
    # Copy & select relevant columns
    bars = (
        df.copy()[["ranks", "analysis_m", "forecast_m", "ens_init_m", "io_m","io_i_m", "io_f_m", "io_a_m"]]
        .sort_values("ranks")
    )

    # subtract off IO from analysis/forecast/ens_init if detailed IO columns are available
    # if not bars[["io_i_m", "io_f_m", "io_a_m"]].isnull().all().all():
    #     bars["analysis_m"] = bars["analysis_m"] - bars["io_a_m"].fillna(0.0)
    #     bars["forecast_m"] = bars["forecast_m"] - bars["io_f_m"].fillna(0.0)
    #     bars["ens_init_m"] = bars["ens_init_m"] - bars["io_i_m"].fillna(0.0)
    #     bars.loc[bars["analysis_m"] < 0.0, "analysis_m"] = 0.0
    #     bars.loc[bars["forecast_m"] < 0.0, "forecast_m"] = 0.0
    #     bars.loc[bars["ens_init_m"] < 0.0, "ens_init_m"] = 0.0

    # Total ranks for labeling
    ranks = bars["ranks"].to_numpy()
    ranks_total = ranks * (model_np + 1)
    x = np.arange(len(ranks_total))  # equal spacing

    width = 0.65

    # Compute totals and percentages
    totals = bars.sum(axis=1).to_numpy(dtype=float)
    bars_pct = bars.div(totals, axis=0) * 100  # each row = 100%

    label_map = {
        "analysis_m": "Analysis",
        "forecast_m": "Forecast",
        "ens_init_m": "Ens Init",
        "io_m": "I/O",
    }

    fig, ax = plt.subplots(figsize=(12, 6.5), dpi=160)

    bottom = np.zeros_like(x, dtype=float)
    for col in ["analysis_m", "forecast_m", "ens_init_m", "io_m"]:
        vals = bars_pct[col].to_numpy(float)
        ax.bar(
            x, vals, width=width, bottom=bottom,
            label=label_map[col], edgecolor="black", linewidth=0.6, alpha=0.9
        )
        bottom += vals

    # Styling
    ax.set_xlabel("Ranks (MPI processes)", fontdict={"fontsize": 18, "fontweight": "bold"})
    # ax.set_ylabel("Time Breakdown (%)",   fontdict={"fontsize": 18, "fontweight": "bold"})
    ax.set_ylabel("Time (mins)",   fontdict={"fontsize": 18, "fontweight": "bold"})
    ax.set_title(title,                   fontdict={"fontsize": 16, "fontweight": "bold"})
    ax.grid(True, axis="y", linestyle="--", alpha=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(ranks_total, fontsize=18, fontweight="bold")
    ax.tick_params(axis="y", labelsize=18)
    ax.set_yticklabels([f"{N:.0f}" for N in ax.get_yticks()], fontsize=18, fontweight="bold")
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')

    # overwrite y-ticks to show minutes instead of percentage
    yticks = ax.get_yticks()
    yticklabels = [f"{(N/100)*totals.max():.0f}" for N in yticks]
    ax.set_yticklabels(yticklabels, fontsize=18, fontweight="bold")
    
    # ax.set_ylim(0, 105)
    ax.legend(loc="upper right", prop={"size": 18, "weight": "bold"})
    ax.margins(x=0.06)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=200)
    plt.close(fig)

def plot_stacked_bars(df: pd.DataFrame, title: str, out_path: Path, model_np=1,flag=None):
    # Copy & select columns
    bars = (
        df.copy()[["ranks", "analysis_m", "forecast_m", "ens_init_m", "assim_m", "io_m"
                   , "io_a_m", "io_f_m", "io_i_m"]]
        .sort_values("ranks")
    )

    # take off I/O from analysis, forecast, ens_init
    bars["analysis_m"] = bars["analysis_m"] - bars["io_a_m"] 
    bars["forecast_m"] = bars["forecast_m"] - bars["io_f_m"]
    bars["ens_init_m"] = bars["ens_init_m"] - bars["io_i_m"]
    bars['assim_m'] = bars['assim_m'] - bars['io_m'] 
   
    # bars = bars[["ranks", "analysis_m", "forecast_m", "ens_init_m", "io_m"]]
    # # bars = bars[bars["analysis_m"] >= 0]
    # # bars = bars[bars["forecast_m"] >= 0]
    # # bars = bars[bars["ens_init_m"] >= 0]
    # # bars = bars[bars["io_m"] >= 0]

    # Use evenly spaced x positions; show actual ranks as labels
    ranks = bars["ranks"].to_numpy()
    ranks_total = ranks * (model_np + 1)      # label to show
    x = np.arange(len(ranks_total))           # 0,1,2,... ensures equal spacing
    width = 0.65                               # bar width in "index" units

    fig, ax = plt.subplots(figsize=(11, 6.5), dpi=160)

    bottom = np.zeros_like(x, dtype=float)
    nice = {
        "forecast_m": "Forecast",
        "analysis_m": "Analysis",
        "ens_init_m": "Ens Init",
        "io_m": "I/O",
    }
    colors = {'forecast_m':'#ff7f0e', 'analysis_m':"#184ec4", 'ens_init_m':'#2ca02c', 'io_m':'#d62728'}
    for col in ["forecast_m", "analysis_m", "ens_init_m", "io_m"]:
        vals = bars[col].to_numpy()
        ax.bar(
            x, vals, width=width, bottom=bottom, label=nice[col], color=colors[col],
            edgecolor="black", linewidth=0.6, alpha=0.9
        )
        bottom += vals

    # Axis cosmetics
    ax.set_xlabel("Ranks (MPI processes)", fontdict={"fontsize": 20, "fontweight": "bold"})
    ax.set_ylabel("Time (minutes)",        fontdict={"fontsize": 20, "fontweight": "bold"})
    # ax.set_title(title,                    fontdict={"fontsize": 20, "fontweight": "bold"})
    ax.grid(True, axis="y", linestyle="--", alpha=0.6)

    # Equal spacing: ticks at the evenly spaced indices, labels are the true ranks
    ax.set_xticks(x)
    ax.set_xticklabels(ranks_total, fontsize=20, fontweight="bold")
    ax.tick_params(axis="y", labelsize=20)
    ax.set_yticklabels([f"{N:.0f}" for N in ax.get_yticks()], fontsize=20, fontweight="bold")

    if flag =='weak_scaling':
        # put the legend outside the plot
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), prop={"size": 18, "weight": "bold"})
    else:   
        ax.legend(loc="upper right", prop={"size": 18, "weight": "bold"})
    ax.margins(x=0.06)                     # padding at plot edges
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

def analyze_icesee(log_path: str, scaling_type: str = "strong_scaling", model_np: int = 1):
    out_dir = Path("icesee_perf_outputs")
    out_dir.mkdir(exist_ok=True, parents=True)

    df = parse_icesee_log(log_path, model_np=model_np)
    if df.empty:
        raise RuntimeError("No timing blocks found in the log. Check formatting or path.")

    # --- Only do wallclock, forecast, analysis, & I/O     speedup/efficiency ---
    wall_metrics = compute_scaling(df, scaling_type, "wall_m")
    fcast_metrics = compute_scaling(df, scaling_type, "forecast_m")
    analysis_metrics = compute_scaling(df, scaling_type, "analysis_m")
    io_metrics = compute_scaling(df, scaling_type, "io_m")

    # Save CSVs
    wall_csv = out_dir / f"{scaling_type}_wallclock_scaling.csv"
    fcast_csv = out_dir / f"{scaling_type}_forecast_scaling.csv"
    analysis_csv = out_dir / f"{scaling_type}_analysis_scaling.csv"
    io_csv = out_dir / f"{scaling_type}_io_scaling.csv"
    df_csv = out_dir / f"{scaling_type}_raw_times.csv"
    wall_metrics.to_csv(wall_csv, index=False)
    fcast_metrics.to_csv(fcast_csv, index=False)
    analysis_metrics.to_csv(analysis_csv, index=False)
    io_metrics.to_csv(io_csv, index=False)
    df.to_csv(df_csv, index=False)

    # --- Plots: speedup & efficiency for each (wallclock and forecast) ---
    # Wallclock
    # Ideal efficiency = 100%; Ideal speedup = N/N0
    N = wall_metrics["ranks"].values
    N0 = wall_metrics["baseline_ranks"].iloc[0] if "baseline_ranks" in wall_metrics else N.min()
    ideal_speedup = N / N0
    plot_line(N, wall_metrics["speedup"].values,
              "Ranks (MPI processes)", "Speedup",
              f"{scaling_type.replace('_',' ').title()}: Wallclock Speedup",
              out_dir / f"{scaling_type}_wallclock_speedup.png",
              ideal=ideal_speedup, logx=False, logy=False, _loglog=True, model_np=model_np)
    ideal_eff = [100.0]*len(N)
    plot_line(N, wall_metrics["efficiency"].values,
              "Ranks (MPI processes)", "Efficiency (%)",
              f"{scaling_type.replace('_',' ').title()}: Wallclock Efficiency",
              out_dir / f"{scaling_type}_wallclock_efficiency.png",
              ideal=ideal_eff, logx=True, logy=False, model_np=model_np)

    # Forecast
    Nf = fcast_metrics["ranks"].values
    N0f = fcast_metrics["baseline_ranks"].iloc[0] if "baseline_ranks" in fcast_metrics else Nf.min()
    ideal_speedup_f = Nf / N0f
    plot_line(Nf, fcast_metrics["speedup"].values,
              "Ranks (MPI processes)", "Speedup",
              f"{scaling_type.replace('_',' ').title()}: Forecast Speedup",
              out_dir / f"{scaling_type}_forecast_speedup.png",
              ideal=ideal_speedup_f, logx=False, logy=False, _loglog=True, model_np=model_np)
    ideal_eff_f = [100.0]*len(Nf)
    plot_line(Nf, fcast_metrics["efficiency"].values,
              "Ranks (MPI processes)", "Efficiency (%)",
              f"{scaling_type.replace('_',' ').title()}: Forecast Efficiency",
              out_dir / f"{scaling_type}_forecast_efficiency.png",
              ideal=ideal_eff_f, logx=True, logy=False, model_np=model_np)
    
    # analysis
    N0a = analysis_metrics["baseline_ranks"].iloc[0] if "baseline_ranks" in analysis_metrics else Nf.min()
    ideal_speedup_a = Nf / N0a
    plot_line(Nf, analysis_metrics["speedup"].values,
              "Ranks (MPI processes)", "Speedup",
              f"{scaling_type.replace('_',' ').title()}: Analysis Speedup",
              out_dir / f"{scaling_type}_analysis_speedup.png",
              ideal=ideal_speedup_a, logx=False, logy=False, _loglog=True, model_np=model_np)
    ideal_eff_a = [100.0]*len(Nf)
    plot_line(Nf, analysis_metrics["efficiency"].values,
              "Ranks (MPI processes)", "Efficiency (%)",
              f"{scaling_type.replace('_',' ').title()}: Analysis Efficiency",
              out_dir / f"{scaling_type}_analysis_efficiency.png",
              ideal=ideal_eff_a, logx=True, logy=False, model_np=model_np)
    
    # I/O
    N0i = io_metrics["baseline_ranks"].iloc[0] if "baseline_ranks" in io_metrics else Nf.min()
    ideal_speedup_i = Nf / N0i
    plot_line(Nf, io_metrics["speedup"].values,
              "Ranks (MPI processes)", "Speedup",
              f"{scaling_type.replace('_',' ').title()}: I/O Speedup",
              out_dir / f"{scaling_type}_io_speedup.png",
              ideal=ideal_speedup_i, logx=False, logy=False, _loglog=True, model_np=model_np)
    ideal_eff_i = [100.0]*len(Nf)
    plot_line(Nf, io_metrics["efficiency"].values,
              "Ranks (MPI processes)", "Efficiency (%)",
              f"{scaling_type.replace('_',' ').title()}: I/O Efficiency",
              out_dir / f"{scaling_type}_io_efficiency.png",
              ideal=ideal_eff_i, logx=True, logy=False, model_np=model_np)

    # --- Stacked bars ---
    plot_stacked_bars(df, f"{scaling_type.replace('_', ' ').title()}: Time Breakdown", out_dir / f"{scaling_type}_stacked_breakdown.png", model_np=model_np,flag=scaling_type)

    return {
        "outputs_dir": str(out_dir.resolve()),
        "csvs": [str(wall_csv), str(fcast_csv), str(df_csv)],
        "plots": [str(p) for p in sorted(out_dir.glob("*.png"))]
    }

# Run on your uploaded weak log as an example (weak scaling).
# manifest = analyze_icesee("weak_timing.log", scaling_type="weak_scaling", model_np=1)
# manifest