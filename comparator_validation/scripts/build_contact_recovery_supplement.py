#!/usr/bin/env python3
"""Build Supplementary Figure S3 from corrected case-study recovery tables."""
from __future__ import annotations

import argparse
import math
from io import BytesIO
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = ROOT / "outputs" / "case_study"
DEFAULT_ENV_PATH = DEFAULT_OUT_DIR / "case_study_environment_recovery.tsv"
DEFAULT_BASIC_PATH = DEFAULT_OUT_DIR / "case_study_crystal_basic_contact_recovery.tsv"

WIDTH, HEIGHT = 3600, 2700
DPI = 300
VARIANTS: Sequence[str] = [
    "PhosphoFill top1",
    "PhosphoFill best-of-3",
    "PyTMs default",
    "PyTMs optimized",
    "PTM-Psi",
]
DISPLAY_LABELS = {
    "PhosphoFill top1": "PhosphoFill top-1",
    "PhosphoFill best-of-3": "PhosphoFill best-of-3",
    "PyTMs default": "PyTMs default",
    "PyTMs optimized": "PyTMs optimised",
    "PTM-Psi": "PTM-Psi",
}

TEXT = "#151515"
MUTED = "#555555"
BORDER = "#d2d2d2"
GRID = "#e7e7e7"
HEADER_BG = "#f4f4f4"
ROW_ALT = "#fafafa"
GREEN = "#cfead1"
YELLOW = "#fff1bd"
RED = "#f4c7c3"
NA_FILL = "#eeeeee"


def font(size: int, bold: bool = False):
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def draw_text(
    draw: ImageDraw.ImageDraw,
    xy: Tuple[float, float],
    value,
    size: int = 32,
    fill: str = TEXT,
    anchor: str = "mm",
    bold: bool = False,
):
    draw.text(xy, str(value), font=font(size, bold), fill=fill, anchor=anchor)


def safe_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def pct_label(recovered, total, recall=None) -> str:
    recovered_i = int(recovered) if not pd.isna(recovered) else 0
    total_i = int(total) if not pd.isna(total) else 0
    if total_i <= 0:
        return "NA"
    if recall is None or math.isnan(safe_float(recall)):
        pct = int(round(100 * recovered_i / total_i))
    else:
        pct = int(round(100 * safe_float(recall)))
    return f"{recovered_i}/{total_i} ({pct}%)"


def recovery_fill(recovered, total, recall=None) -> str:
    total_f = safe_float(total)
    if math.isnan(total_f) or total_f <= 0:
        return NA_FILL
    if recall is None or math.isnan(safe_float(recall)):
        ratio = safe_float(recovered) / total_f
    else:
        ratio = safe_float(recall)
    if ratio >= 0.75:
        return GREEN
    if ratio >= 0.4:
        return YELLOW
    return RED


def compact_case_label(case_id: str) -> str:
    parts = str(case_id).split()
    if len(parts) >= 4:
        return f"{' '.join(parts[:3])} ({parts[3]})"
    return str(case_id)


def draw_case_label(draw: ImageDraw.ImageDraw, x: int, y: int, case_id: str):
    draw_text(draw, (x, y), compact_case_label(case_id), size=28, anchor="lm")


def panel_rows(df: pd.DataFrame) -> List[str]:
    return list(dict.fromkeys(df["case_id"]))


def panel_frame(draw: ImageDraw.ImageDraw, box: Tuple[int, int, int, int], label: str, title: str, subtitle: str):
    x, y, w, h = box
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, outline=BORDER, width=4, fill="white")
    draw_text(draw, (x + 34, y + 40), label, size=58, bold=True, anchor="lm")
    draw_text(draw, (x + w / 2, y + 42), title, size=44, bold=True)
    draw_text(draw, (x + w / 2, y + 88), subtitle, size=25, fill=MUTED)


def draw_heatmap_panel(
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    label: str,
    title: str,
    subtitle: str,
    df: pd.DataFrame,
    total_col: str,
    recovered_col: str,
    recall_col: str,
    total_header: str,
):
    panel_frame(draw, box, label, title, subtitle)
    x, y, w, h = box
    rows = panel_rows(df)
    left = x + 46
    top = y + 175
    case_w = 375
    total_w = 145
    gap = 18
    table_w = w - 92
    cell_w = (table_w - case_w - total_w - gap * (len(VARIANTS) + 1)) / len(VARIANTS)
    row_h = 65

    draw.rectangle((left, top - 44, left + table_w, top + 18), fill=HEADER_BG)
    draw_text(draw, (left, top - 14), "Case", size=27, anchor="lm", bold=True)
    draw_text(draw, (left + case_w + total_w / 2, top - 14), total_header, size=27, bold=True)
    variant_x0 = left + case_w + total_w + gap
    for j, variant in enumerate(VARIANTS):
        cx = variant_x0 + j * (cell_w + gap) + cell_w / 2
        draw_text(draw, (cx, top - 14), DISPLAY_LABELS.get(variant, variant), size=25, bold=True)

    for i, case_id in enumerate(rows):
        yy = top + 44 + i * row_h
        sub = df[df["case_id"] == case_id]
        first = sub.iloc[0]
        if i % 2:
            draw.rectangle((left, yy - 31, left + table_w, yy + 31), fill=ROW_ALT)
        draw_case_label(draw, left, yy, case_id)
        total = int(first[total_col])
        draw_text(draw, (left + case_w + total_w / 2, yy), str(total), size=29, bold=True)
        for j, variant in enumerate(VARIANTS):
            row = sub[sub["variant"] == variant]
            cell_x = variant_x0 + j * (cell_w + gap)
            fill = NA_FILL
            label_text = "NA"
            if not row.empty:
                r = row.iloc[0]
                fill = recovery_fill(r[recovered_col], r[total_col], r[recall_col])
                label_text = pct_label(r[recovered_col], r[total_col], r[recall_col])
            draw.rounded_rectangle(
                (cell_x, yy - 27, cell_x + cell_w, yy + 27),
                radius=8,
                fill=fill,
                outline="#dddddd",
                width=2,
            )
            draw_text(draw, (cell_x + cell_w / 2, yy), label_text, size=27)


def legend(draw: ImageDraw.ImageDraw, x: int, y: int):
    draw_text(draw, (x, y), "Cell color:", size=27, bold=True, anchor="lm")
    items = [
        (GREEN, ">=75% recovered"),
        (YELLOW, "40-74% recovered"),
        (RED, "<40% recovered"),
    ]
    lx = x + 155
    for fill, label in items:
        draw.rounded_rectangle((lx, y - 20, lx + 54, y + 20), radius=7, fill=fill, outline="#dddddd")
        draw_text(draw, (lx + 70, y), label, size=25, anchor="lm", fill=MUTED)
        lx += 335
    draw_text(
        draw,
        (x + 1300, y),
        "All cells report recovered crystal contacts / total crystal contacts.",
        size=25,
        anchor="lm",
        fill=MUTED,
    )


def build_summary_figure(env: pd.DataFrame, basic: pd.DataFrame) -> Image.Image:
    image = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (WIDTH / 2, 62), "Supplementary contact-recovery summary", size=54, bold=True)
    draw_text(
        draw,
        (WIDTH / 2, 112),
        "Compact recovery heatmaps for whole local environments, basic residues, and donor atom-phosphate oxygen contacts.",
        size=29,
        fill=MUTED,
    )

    margin_x = 70
    panel_w = WIDTH - 2 * margin_x
    panel_h = 715
    draw_heatmap_panel(
        draw,
        (margin_x, 160, panel_w, panel_h),
        "(a)",
        "5 A crystal environment contact recovery",
        "All non-water heavy-atom contacts around the phosphosite marker atoms.",
        env,
        "phospho_crystal_contact_count_5A",
        "recovered_crystal_contact_count_5A",
        "contact_recall_5A",
        "Crystal contacts",
    )
    draw_heatmap_panel(
        draw,
        (margin_x, 905, panel_w, panel_h),
        "(b)",
        "4 A basic-residue contact recovery",
        "Recovery of crystal Arg/Lys/His residues contacting phosphate oxygens.",
        basic,
        "phospho_crystal_basic_contact_count_4A",
        "recovered_crystal_basic_contact_count_4A",
        "basic_contact_recall_4A",
        "Crystal basic",
    )
    draw_heatmap_panel(
        draw,
        (margin_x, 1650, panel_w, panel_h),
        "(c)",
        "4 A donor atom-phosphate oxygen recovery",
        "Recovery of donor atom contacts to any phosphate oxygen.",
        basic,
        "phospho_crystal_basic_atomic_contact_count_4A",
        "recovered_crystal_basic_atomic_contact_count_4A",
        "basic_atomic_contact_recall_4A",
        "Crystal donors",
    )
    legend(draw, margin_x + 30, 2440)
    draw_text(
        draw,
        (margin_x + 30, 2502),
        "Detailed residue and donor-atom labels are provided in figs3_contact_recovery_detail.tsv.",
        size=25,
        anchor="lm",
        fill=MUTED,
    )
    return image


DETAIL_COLUMNS = [
    "protein",
    "site_label",
    "case_id",
    "pdb_id",
    "tool",
    "variant",
    "phospho_crystal_contact_count_5A",
    "predicted_contact_count_5A",
    "recovered_crystal_contact_count_5A",
    "missing_crystal_contact_count_5A",
    "extra_predicted_contact_count_5A",
    "contact_recall_5A",
    "contact_precision_5A",
    "phospho_crystal_contact_set_5A",
    "predicted_contact_set_5A",
    "recovered_crystal_contacts_5A",
    "missing_crystal_contacts_5A",
    "extra_predicted_contacts_5A",
    "phospho_crystal_basic_contact_count_4A",
    "predicted_basic_contact_count_4A",
    "recovered_crystal_basic_contact_count_4A",
    "missing_crystal_basic_contact_count_4A",
    "extra_predicted_basic_contact_count_4A",
    "basic_contact_recall_4A",
    "basic_contact_precision_4A",
    "phospho_crystal_basic_contact_set_4A",
    "phospho_crystal_basic_contact_set_4A_mapped",
    "predicted_basic_contact_set_4A",
    "recovered_crystal_basic_contacts_4A",
    "missing_crystal_basic_contacts_4A",
    "extra_predicted_basic_contacts_4A",
    "phospho_crystal_basic_atomic_contact_count_4A",
    "predicted_basic_atomic_contact_count_4A",
    "basic_atomic_contact_delta_4A",
    "recovered_crystal_basic_atomic_contact_count_4A",
    "missing_crystal_basic_atomic_contact_count_4A",
    "extra_predicted_basic_atomic_contact_count_4A",
    "basic_atomic_contact_recall_4A",
    "basic_atomic_contact_precision_4A",
    "phospho_crystal_basic_atomic_contact_set_4A",
    "phospho_crystal_basic_atomic_contact_set_4A_mapped",
    "predicted_basic_atomic_contact_set_4A",
    "recovered_crystal_basic_atomic_contacts_4A",
    "missing_crystal_basic_atomic_contacts_4A",
    "extra_predicted_basic_atomic_contacts_4A",
    "structure_path",
]


def build_detail_table(env: pd.DataFrame, basic: pd.DataFrame) -> pd.DataFrame:
    basic_cols = [
        "case_id",
        "variant",
        "phospho_crystal_basic_contact_set_4A_mapped",
        "phospho_crystal_basic_atomic_contact_set_4A_mapped",
    ]
    merged = env.merge(
        basic[basic_cols],
        on=["case_id", "variant"],
        how="left",
        suffixes=("", "_from_basic"),
    )
    for col in [
        "phospho_crystal_basic_contact_set_4A_mapped",
        "phospho_crystal_basic_atomic_contact_set_4A_mapped",
    ]:
        alt = f"{col}_from_basic"
        if alt in merged:
            merged[col] = merged[col].fillna(merged[alt])
    available = [col for col in DETAIL_COLUMNS if col in merged.columns]
    table = merged[available].copy()
    variant_order = {variant: i for i, variant in enumerate(VARIANTS)}
    table["_variant_order"] = table["variant"].map(variant_order).fillna(99)
    table.sort_values(["protein", "position" if "position" in table.columns else "case_id", "pdb_id", "_variant_order"], inplace=True)
    table.drop(columns=["_variant_order"], inplace=True)
    return table


def build_summary_long(env: pd.DataFrame, basic: pd.DataFrame) -> pd.DataFrame:
    panels = [
        ("a", "5A_environment", env, "phospho_crystal_contact_count_5A", "recovered_crystal_contact_count_5A", "contact_recall_5A"),
        ("b", "4A_basic_residue", basic, "phospho_crystal_basic_contact_count_4A", "recovered_crystal_basic_contact_count_4A", "basic_contact_recall_4A"),
        ("c", "4A_donor_atom", basic, "phospho_crystal_basic_atomic_contact_count_4A", "recovered_crystal_basic_atomic_contact_count_4A", "basic_atomic_contact_recall_4A"),
    ]
    rows: List[Dict[str, object]] = []
    for panel_id, metric, df, total_col, recovered_col, recall_col in panels:
        for _, row in df.iterrows():
            total = int(row[total_col])
            recovered = int(row[recovered_col])
            recall = safe_float(row[recall_col])
            rows.append({
                "panel": panel_id,
                "metric": metric,
                "protein": row.get("protein", ""),
                "site_label": row.get("site_label", ""),
                "case_id": row["case_id"],
                "pdb_id": row.get("pdb_id", ""),
                "tool": row.get("tool", ""),
                "variant": row["variant"],
                "recovered": recovered,
                "total": total,
                "percent_recovered": 100 * (recall if not math.isnan(recall) else recovered / total) if total else float("nan"),
            })
    return pd.DataFrame(rows)


def save_pdf(image: Image.Image, path: Path):
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfgen import canvas

    page_width_pt = WIDTH / DPI * 72
    page_height_pt = HEIGHT / DPI * 72
    buffer = BytesIO()
    image.save(buffer, format="PNG", dpi=(DPI, DPI))
    buffer.seek(0)
    pdf = canvas.Canvas(str(path), pagesize=(page_width_pt, page_height_pt))
    pdf.drawImage(ImageReader(buffer), 0, 0, width=page_width_pt, height=page_height_pt)
    pdf.showPage()
    pdf.save()


def assert_inputs_exist(paths: Iterable[Path]):
    missing = [path for path in paths if not path.exists()]
    if missing:
        joined = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing required input files:\n{joined}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env", type=Path, default=DEFAULT_ENV_PATH, help="Corrected environment-recovery TSV")
    parser.add_argument("--basic", type=Path, default=DEFAULT_BASIC_PATH, help="Corrected basic-contact-recovery TSV")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    return parser.parse_args()


def write_caption(path: Path) -> None:
    text = r"""\caption*{\textbf{Supplementary Figure S3.} Contact-recovery summary across CDK2 and ERK2 phosphosite case studies. Heatmaps compare PhosphoFill top-1, PhosphoFill best-of-3, PyTMs default, PyTMs optimised, and PTM-Psi for CDK2 pThr160 and ERK2 pThr185/pTyr187 strip-and-regraft cases. PhosphoFill best-of-3 denotes the independently minimised rank with the lowest phosphate RMSD to the corresponding experimental coordinates. Each cell reports recovered crystal contacts / total crystal contacts; green indicates $\geq75\%$ recovery, yellow 40--74\%, and red $<40\%$. \textbf{(a)} Recovery of local non-water heavy-atom contacts within 5~\AA\ of the phosphosite marker atoms. PhosphoFill top-1 recovers all or nearly all crystal-environment contacts across the six cases, whereas comparator tools show partial losses, especially for ERK2 pTyr187. \textbf{(b)} Recovery of basic Arg/Lys/His residue contacts within 4~\AA\ of phosphate oxygens. PhosphoFill top-1 and best-of-3 recover 100\% of crystal basic-residue contacts in all six cases; comparator methods show case-dependent losses. \textbf{(c)} Recovery of donor atom--phosphate oxygen contacts within 4~\AA. PhosphoFill preserves all crystal donor-atom contacts across all six cases, while comparator methods lose one or more donor contacts in multiple examples, including complete loss in some ERK2 pTyr187 predictions.}
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    assert_inputs_exist([args.env, args.basic])
    env = pd.read_csv(args.env, sep="\t")
    basic = pd.read_csv(args.basic, sep="\t")
    env = env[env["variant"].isin(VARIANTS)].copy()
    basic = basic[basic["variant"].isin(VARIANTS)].copy()

    figure_dir = args.outdir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    output_png = figure_dir / "figs3_contact_recovery.png"
    output_tif = figure_dir / "figs3_contact_recovery.tif"
    output_pdf = figure_dir / "figs3_contact_recovery.pdf"
    output_table = args.outdir / "figs3_contact_recovery_detail.tsv"
    output_source = args.outdir / "figs3_contact_recovery_source_values.tsv"

    image = build_summary_figure(env, basic)
    image.save(output_png, dpi=(DPI, DPI))
    image.save(output_tif, dpi=(DPI, DPI), compression="tiff_lzw")
    save_pdf(image, output_pdf)

    detail = build_detail_table(env, basic)
    detail.to_csv(output_table, sep="\t", index=False)

    summary_long = build_summary_long(env, basic)
    summary_long["display_variant"] = summary_long["variant"].map(DISPLAY_LABELS)
    summary_long.to_csv(output_source, sep="\t", index=False)
    write_caption(args.outdir / "figs3_contact_recovery_caption.tex")

    print(output_png)
    print(output_tif)
    print(output_pdf)
    print(output_table)
    print(image.size)
    print(f"detail_rows={len(detail)} detail_cols={len(detail.columns)}")


if __name__ == "__main__":
    main()
