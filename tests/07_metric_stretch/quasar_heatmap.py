import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# PROJECT COEUS: QUASAR MEAN-REDSHIFT SECTOR AUDIT
# Fixes:
#   - numeric bin centers, not categorical x-axis
#   - per-sector counts printed/saved
#   - optional minimum count threshold to avoid noisy bins
# ============================================================

DATA_DIR = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"
QSO_PATH = os.path.join(DATA_DIR, "dr16q_lean_audit.csv")
OUTPUT_PNG = os.path.join(DATA_DIR, "quasar_mean_z_audit_fixed.png")
OUTPUT_CSV = os.path.join(DATA_DIR, "quasar_mean_z_audit_fixed_summary.csv")

Z_MIN = 1.0
Z_MAX = 2.5
BIN_WIDTH_DEG = 15
TARGET_CENTER = 240
TARGET_HALF_WIDTH = 20
MIN_COUNT_PER_BIN = 10


def main() -> None:
    if not os.path.exists(QSO_PATH):
        raise FileNotFoundError(QSO_PATH)

    df = pd.read_csv(QSO_PATH)
    df.columns = [c.upper().strip() for c in df.columns]
    if not {"RA", "Z"}.issubset(df.columns):
        raise ValueError("Expected columns RA and Z")

    df = df.dropna(subset=["RA", "Z"]).copy()
    df["RA"] = df["RA"] % 360.0
    df = df[(df["Z"] >= Z_MIN) & (df["Z"] <= Z_MAX)].copy()

    # Build explicit numeric sectors
    starts = np.arange(0, 360, BIN_WIDTH_DEG, dtype=float)
    centers = (starts + BIN_WIDTH_DEG / 2.0) % 360.0

    means = []
    counts = []
    labels = []
    for start, center in zip(starts, centers):
        end = (start + BIN_WIDTH_DEG) % 360.0
        if start < end:
            mask = (df["RA"] >= start) & (df["RA"] < end)
        else:
            mask = (df["RA"] >= start) | (df["RA"] < end)
        sub = df.loc[mask, "Z"]
        n = int(sub.shape[0])
        counts.append(n)
        means.append(float(sub.mean()) if n >= MIN_COUNT_PER_BIN else np.nan)
        labels.append(f"{int(start)}-{int((start + BIN_WIDTH_DEG) % 360)}")

    out = pd.DataFrame({
        "ra_bin_start_deg": starts,
        "ra_bin_center_deg": centers,
        "ra_bin_label": labels,
        "count": counts,
        "mean_z": means,
    })
    out.to_csv(OUTPUT_CSV, index=False)

    fig, ax1 = plt.subplots(figsize=(14, 6))
    ax1.plot(centers, means, marker="o", linestyle="-", color="purple", linewidth=2, label="Mean quasar redshift")
    ax1.axvline(TARGET_CENTER, color="red", linestyle="--", label=f"The Basin ({TARGET_CENTER}°)")
    ax1.axvspan(TARGET_CENTER - TARGET_HALF_WIDTH, TARGET_CENTER + TARGET_HALF_WIDTH,
                color="red", alpha=0.1, label="Target sector")
    ax1.set_title("METRIC STRETCH AUDIT: Redshift Consistency Check", fontsize=16)
    ax1.set_xlabel("Right Ascension Sector Center (deg)", fontsize=12)
    ax1.set_ylabel("Mean Redshift (Z) per Sector", fontsize=12)
    ax1.set_xlim(0, 360)
    ax1.set_xticks(np.arange(0, 361, 30))
    ax1.grid(alpha=0.3)

    # Secondary axis for counts so sparse bins are obvious
    ax2 = ax1.twinx()
    ax2.bar(centers, counts, width=BIN_WIDTH_DEG * 0.8, color="gray", alpha=0.15, label="Counts per sector")
    ax2.set_ylabel("Quasar Count per Sector", fontsize=12)

    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper right")

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
    plt.show()

    print(f"Saved plot: {OUTPUT_PNG}")
    print(f"Saved summary: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
