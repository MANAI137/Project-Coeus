import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# PROJECT COEUS: QUASAR MANIFOLD AUDIT (FOOTPRINT-PROXY NORMALIZED)
# Methodology:
#   - No random catalogs
#   - No isotropic priors
#   - Denominator is built from the full observed quasar catalog footprint,
#     not from the shell being tested
# ============================================================

DATA_DIR = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"
QSO_PATH = os.path.join(DATA_DIR, "dr16q_lean_audit.csv")
OUTPUT_PNG = os.path.join(DATA_DIR, "quasar_manifold_audit_fixed.png")
OUTPUT_CSV = os.path.join(DATA_DIR, "quasar_manifold_audit_fixed_summary.csv")

SHELLS = [
    {"z_min": 0.1, "z_max": 0.8, "title": "Q-SHELL 1: Local Overlap",      "color": "cyan"},
    {"z_min": 0.8, "z_max": 1.5, "title": "Q-SHELL 2: The Acceleration",   "color": "blue"},
    {"z_min": 1.5, "z_max": 2.2, "title": "Q-SHELL 3: The Deep Bridge",    "color": "darkblue"},
    {"z_min": 2.2, "z_max": 4.0, "title": "Q-SHELL 4: The Terminal Horizon","color": "purple"},
]

SECTOR_CENTERS = np.array([15, 60, 105, 150, 195, 240, 285, 330], dtype=float)
SECTOR_HALF_WIDTH = 22.5
PIXEL_SIZE_DEG = 1.0
HIGHLIGHT_SECTORS = {195, 240}
MIN_AREA_SR = 0.05


def standardize(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [c.upper().strip() for c in out.columns]
    required = {"RA", "DEC", "Z"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    out["RA"] = out["RA"] % 360.0
    return out.dropna(subset=["RA", "DEC", "Z"]).copy()


def in_sector(ra_series: pd.Series, center_deg: float, half_width_deg: float) -> pd.Series:
    delta = ((ra_series - center_deg + 180.0) % 360.0) - 180.0
    return np.abs(delta) < half_width_deg


def build_bulk_pixels(bulk_df: pd.DataFrame) -> pd.DataFrame:
    px = bulk_df.copy()
    px["RA_PIX"] = np.floor(px["RA"]).astype(int) % 360
    px["DEC_PIX"] = np.floor(px["DEC"]).astype(int)
    px = px.drop_duplicates(subset=["RA_PIX", "DEC_PIX"]).copy()
    return px


def pixel_solid_angle_sr(dec_pix_deg: pd.Series, pixel_size_deg: float = 1.0) -> np.ndarray:
    half = pixel_size_deg / 2.0
    dec_center = dec_pix_deg.astype(float) + 0.5
    dec_lo = np.deg2rad(dec_center - half)
    dec_hi = np.deg2rad(dec_center + half)
    d_ra = np.deg2rad(pixel_size_deg)
    return d_ra * (np.sin(dec_hi) - np.sin(dec_lo))


def footprint_area_from_bulk(bulk_df: pd.DataFrame, center_deg: float) -> float:
    px = build_bulk_pixels(bulk_df)
    px = px[in_sector(px["RA_PIX"] + 0.5, center_deg, SECTOR_HALF_WIDTH)].copy()
    if px.empty:
        return np.nan
    return float(pixel_solid_angle_sr(px["DEC_PIX"], PIXEL_SIZE_DEG).sum())


def shell_density(shell_df: pd.DataFrame, bulk_df: pd.DataFrame, center_deg: float) -> float:
    area = footprint_area_from_bulk(bulk_df, center_deg)
    if not np.isfinite(area) or area <= MIN_AREA_SR:
        return np.nan
    n_shell = int(in_sector(shell_df["RA"], center_deg, SECTOR_HALF_WIDTH).sum())
    return n_shell / area


def main() -> None:
    if not os.path.exists(QSO_PATH):
        raise FileNotFoundError(QSO_PATH)

    qso_all = standardize(pd.read_csv(QSO_PATH))

    results = []
    fig = plt.figure(figsize=(20, 18))
    angles = np.deg2rad(SECTOR_CENTERS)
    width = np.deg2rad(45.0)

    for i, shell in enumerate(SHELLS, start=1):
        q_shell = qso_all[(qso_all["Z"] >= shell["z_min"]) & (qso_all["Z"] <= shell["z_max"])].copy()

        # Footprint sync: use the full quasar catalog, but only inside the DEC strip
        # occupied by the shell. This keeps denominator and numerator on the same
        # latitudinal footprint without using the shell itself as the area source.
        dec_min = q_shell["DEC"].min()
        dec_max = q_shell["DEC"].max()
        q_bulk = qso_all[(qso_all["DEC"] >= dec_min) & (qso_all["DEC"] <= dec_max)].copy()

        row_out = []
        densities = []
        for center in SECTOR_CENTERS:
            rho = shell_density(q_shell, q_bulk, center)
            n = int(in_sector(q_shell["RA"], center, SECTOR_HALF_WIDTH).sum())
            densities.append(rho)
            row_out.append({
                "shell_title": shell["title"],
                "z_min": shell["z_min"],
                "z_max": shell["z_max"],
                "sector_center_deg": center,
                "count": n,
                "rho_count_per_sr": rho,
                "shell_total": len(q_shell),
            })
        results.extend(row_out)

        ax = fig.add_subplot(2, 2, i, projection="polar")
        bars = ax.bar(
            angles,
            np.nan_to_num(np.array(densities, dtype=float), nan=0.0),
            width=width,
            color=shell["color"],
            alpha=0.7,
            edgecolor="black",
        )

        for j, center in enumerate(SECTOR_CENTERS):
            if center in HIGHLIGHT_SECTORS:
                bars[j].set_edgecolor("red")
                bars[j].set_linewidth(2.0)

        ymax = np.nanmax(densities) if np.isfinite(np.nanmax(densities)) else 1.0
        for angle, center, rho in zip(angles, SECTOR_CENTERS, densities):
            if np.isfinite(rho):
                ax.text(
                    angle,
                    rho + ymax * 0.08,
                    f"ρ={rho:.0f}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                )

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(angles)
        ax.set_xticklabels([f"{int(c)}°" for c in SECTOR_CENTERS])
        ax.set_title(f"{shell['title']}\nz=[{shell['z_min']}-{shell['z_max']}]", fontsize=15, pad=30)

    summary = pd.DataFrame(results)
    summary.to_csv(OUTPUT_CSV, index=False)

    plt.suptitle(
        "PROJECT COEUS: QUASAR MANIFOLD AUDIT\n"
        "Multi-Tracer Footprint-Proxy Density (No Random Catalogs)",
        fontsize=22,
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
    plt.show()

    print(f"Saved plot: {OUTPUT_PNG}")
    print(f"Saved summary: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
