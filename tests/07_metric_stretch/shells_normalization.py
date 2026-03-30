import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# PROJECT COEUS: BULK-FOOTPRINT NORMALIZED SHELL AUDIT
# No random catalogs. No isotropic priors.
# Density is measured relative to observed survey footprint only.
# ============================================================

DATA_DIR = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"
OUTPUT_PNG = os.path.join(DATA_DIR, "normalized_shell_compass.png")
OUTPUT_CSV = os.path.join(DATA_DIR, "normalized_shell_compass_summary.csv")

SHELL_MAP = [
    {
        "locked": "shell_1_2mrs_locked_audit.csv",
        "bulk":   "2mrs_360_bulk_audit.csv",
        "title":  "SHELL 1: The Launch (2MRS)",
        "color":  "indigo",
        "sync":   {}  # fill if needed
    },
    {
        "locked": "shell_2_6dfgs_locked_audit.csv",
        "bulk":   "6dfgs_360_bulk_audit.csv",
        "title":  "SHELL 2: The Boundary (6dFGS)",
        "color":  "darkorange",
        "sync":   {}  # fill if needed
    },
    {
        "locked": "shell_3_sdss_locked_audit.csv",
        "bulk":   "sdss_360_bulk_audit.csv",
        "title":  "SHELL 3: The Bridge (SDSS)",
        "color":  "darkgreen",
        "sync":   {}  # fill if needed
    },
    {
        "locked": "shell_4_desi_locked_audit.csv",
        "bulk":   "desi_360_bulk_audit.csv",
        "title":  "SHELL 4: The Basin (DESI)",
        "color":  "darkred",
        "sync":   {
            # Example:
            # "dec_min": 25,
            # "dec_max": 42,
            # "ra_ranges": [(135, 255)]
        }
    },
]

SECTOR_CENTERS = np.array([15, 60, 105, 150, 195, 240, 285, 330], dtype=float)
SECTOR_HALF_WIDTH = 22.5
PIXEL_SIZE_DEG = 1.0

# Sectors to visually emphasize as the candidate manifold path
HIGHLIGHT_SECTORS = {150, 195, 240}


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.upper().strip() for c in df.columns]
    required = {"RA", "DEC"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    df["RA"] = df["RA"] % 360.0
    return df


def apply_sync_cuts(df: pd.DataFrame, sync: dict) -> pd.DataFrame:
    """
    Optional footprint-matching cuts.
    Supported keys:
      - dec_min
      - dec_max
      - ra_ranges : list of (ra_min, ra_max), supports wrap
    """
    out = df.copy()

    if "dec_min" in sync:
        out = out[out["DEC"] >= sync["dec_min"]]
    if "dec_max" in sync:
        out = out[out["DEC"] <= sync["dec_max"]]

    if "ra_ranges" in sync and sync["ra_ranges"]:
        masks = []
        for ra_min, ra_max in sync["ra_ranges"]:
            ra_min %= 360.0
            ra_max %= 360.0
            if ra_min <= ra_max:
                masks.append((out["RA"] >= ra_min) & (out["RA"] < ra_max))
            else:
                masks.append((out["RA"] >= ra_min) | (out["RA"] < ra_max))
        ra_mask = np.logical_or.reduce(masks)
        out = out[ra_mask]

    return out


def infer_dec_match_from_shell(shell_df: pd.DataFrame, bulk_df: pd.DataFrame) -> pd.DataFrame:
    """
    Conservative footprint sync: restrict bulk to the shell DEC span.
    This is not a formal survey mask, but it keeps numerator and denominator
    in the same latitudinal strip.
    """
    dec_min = shell_df["DEC"].min()
    dec_max = shell_df["DEC"].max()
    return bulk_df[(bulk_df["DEC"] >= dec_min) & (bulk_df["DEC"] <= dec_max)].copy()


def in_sector(ra_series: pd.Series, center_deg: float, half_width_deg: float) -> pd.Series:
    """
    Returns True if angle lies within wrapped angular distance of sector center.
    """
    delta = ((ra_series - center_deg + 180.0) % 360.0) - 180.0
    return np.abs(delta) < half_width_deg


def build_bulk_pixels(bulk_df: pd.DataFrame) -> pd.DataFrame:
    """
    Construct a 1° x 1° pixelized footprint proxy from the bulk survey.
    """
    px = bulk_df.copy()
    px["RA_PIX"] = np.floor(px["RA"]).astype(int) % 360
    px["DEC_PIX"] = np.floor(px["DEC"]).astype(int)
    px = px.drop_duplicates(subset=["RA_PIX", "DEC_PIX"]).copy()
    return px


def pixel_solid_angle_sr(dec_pix_deg: pd.Series, pixel_size_deg: float = 1.0) -> np.ndarray:
    """
    Exact spherical area for a pixel with height pixel_size_deg in DEC and
    width pixel_size_deg in RA, centered on DEC pixel center.
    Returned in steradians.
    """
    half = pixel_size_deg / 2.0
    dec_center = dec_pix_deg.astype(float) + 0.5
    dec_lo = np.deg2rad(dec_center - half)
    dec_hi = np.deg2rad(dec_center + half)
    d_ra = np.deg2rad(pixel_size_deg)
    return d_ra * (np.sin(dec_hi) - np.sin(dec_lo))


def sector_area_from_bulk_pixels(bulk_pixels: pd.DataFrame, center_deg: float) -> float:
    sector_px = bulk_pixels[in_sector(bulk_pixels["RA_PIX"] + 0.5, center_deg, SECTOR_HALF_WIDTH)].copy()
    if sector_px.empty:
        return np.nan
    return float(pixel_solid_angle_sr(sector_px["DEC_PIX"], PIXEL_SIZE_DEG).sum())


def sector_count(shell_df: pd.DataFrame, center_deg: float) -> int:
    return int(in_sector(shell_df["RA"], center_deg, SECTOR_HALF_WIDTH).sum())


def compute_shell_densities(shell_df: pd.DataFrame, bulk_df: pd.DataFrame) -> pd.DataFrame:
    bulk_pixels = build_bulk_pixels(bulk_df)

    rows = []
    for center in SECTOR_CENTERS:
        n = sector_count(shell_df, center)
        area_sr = sector_area_from_bulk_pixels(bulk_pixels, center)

        if np.isnan(area_sr) or area_sr <= 0:
            rho = np.nan
        else:
            rho = n / area_sr

        rows.append({
            "sector_center_deg": center,
            "count": n,
            "area_sr": area_sr,
            "rho_count_per_sr": rho
        })

    return pd.DataFrame(rows)


def plot_shells(shell_results: list[tuple[dict, pd.DataFrame]]) -> None:
    fig = plt.figure(figsize=(20, 18))
    angles = np.deg2rad(SECTOR_CENTERS)
    width = np.deg2rad(45.0)

    for i, (meta, result_df) in enumerate(shell_results, start=1):
        ax = fig.add_subplot(2, 2, i, projection="polar")

        heights = result_df["rho_count_per_sr"].to_numpy(dtype=float)
        bars = ax.bar(
            angles,
            np.nan_to_num(heights, nan=0.0),
            width=width,
            color=meta["color"],
            alpha=0.7,
            edgecolor="black"
        )

        # Highlight candidate manifold sectors
        for j, center in enumerate(SECTOR_CENTERS):
            if center in HIGHLIGHT_SECTORS:
                bars[j].set_edgecolor("red")
                bars[j].set_linewidth(2.0)

        ymax = np.nanmax(heights) if np.isfinite(np.nanmax(heights)) else 1.0
        for angle, center, rho in zip(angles, SECTOR_CENTERS, heights):
            if np.isfinite(rho):
                ax.text(
                    angle,
                    rho + ymax * 0.08,
                    f"ρ={rho:.1f}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold"
                )

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(angles)
        ax.set_xticklabels([f"{int(x)}°" for x in SECTOR_CENTERS])
        ax.set_title(f"{meta['title']}\nBULK-FOOTPRINT NORMALIZED DENSITY", fontsize=15, pad=30)

    plt.suptitle(
        "PROJECT COEUS: BULK-FOOTPRINT NORMALIZED SHELL AUDIT\n"
        "Discrete Cross-Sections of the 168° → 240° Manifold",
        fontsize=22,
        y=1.02
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
    plt.show()


def main() -> None:
    all_rows = []
    shell_results = []

    for meta in SHELL_MAP:
        locked_path = os.path.join(DATA_DIR, meta["locked"])
        bulk_path = os.path.join(DATA_DIR, meta["bulk"])

        if not os.path.exists(locked_path):
            print(f"Missing locked file: {locked_path}")
            continue
        if not os.path.exists(bulk_path):
            print(f"Missing bulk file: {bulk_path}")
            continue

        shell_df = standardize_columns(pd.read_csv(locked_path))
        bulk_df = standardize_columns(pd.read_csv(bulk_path))

        # First apply any explicit survey-specific sync cuts
        bulk_df = apply_sync_cuts(bulk_df, meta.get("sync", {}))
        shell_df = apply_sync_cuts(shell_df, meta.get("sync", {}))

        # Then force bulk into the shell's DEC strip
        bulk_df = infer_dec_match_from_shell(shell_df, bulk_df)

        result_df = compute_shell_densities(shell_df, bulk_df)
        result_df["shell_title"] = meta["title"]
        result_df["locked_file"] = meta["locked"]
        result_df["bulk_file"] = meta["bulk"]

        all_rows.append(result_df)
        shell_results.append((meta, result_df))

    if not all_rows:
        raise RuntimeError("No valid shell/bulk pairs were processed.")

    summary_df = pd.concat(all_rows, ignore_index=True)
    summary_df.to_csv(OUTPUT_CSV, index=False)
    plot_shells(shell_results)

    print(f"Saved plot: {OUTPUT_PNG}")
    print(f"Saved summary: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()