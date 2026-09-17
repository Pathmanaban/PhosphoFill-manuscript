#!/usr/bin/env python3
"""Build a more readable, non-destructive CDK2 manuscript composite figure."""
from __future__ import annotations

import math
from io import BytesIO
from pathlib import Path
from typing import List, Sequence, Tuple

import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

import case_study_tool_analysis as base


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "case_study_outputs"
FIG_DIR = OUT_DIR / "figures"
OUTPUT_PNG = FIG_DIR / "CDK2_case_study_manuscript_readable_18x12_300dpi.png"
OUTPUT_TIF = FIG_DIR / "CDK2_case_study_manuscript_readable_18x12_300dpi.tif"
OUTPUT_PDF = FIG_DIR / "CDK2_case_study_manuscript_readable_18x12_300dpi.pdf"

WIDTH, HEIGHT = 5400, 3600
DPI = 300
VARIANTS: Sequence[str] = base.SUPPLEMENTARY_VARIANTS
PANEL_COLORS = {
    **base.COLORS,
    "PhosphoFill top1": base.COLORS["PhosphoFill top1"],
    "PhosphoFill best-of-3": base.COLORS["PTM-Psi"],
    "PyTMs default": "#F59127",
    "PyTMs optimized": "#D4A219",
    "PTM-Psi": "#C2185B",
}

PANEL_BG = "#ffffff"
BORDER = "#d0d0d0"
BORDER_WIDTH = 4
TEXT = "#111111"
MUTED = "#555555"
GRID = "#e8e8e8"


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
    size: int = 34,
    fill: str = TEXT,
    anchor: str = "mm",
    bold: bool = False,
):
    draw.text(xy, str(value), font=font(size, bold), fill=fill, anchor=anchor)


def text_width(draw: ImageDraw.ImageDraw, value: str, size: int, bold: bool = False) -> float:
    return draw.textlength(str(value), font=font(size, bold))


def safe_float(value) -> float:
    return base.safe_float(value)


def wrap_case_label(case_id: str) -> List[str]:
    parts = str(case_id).split()
    if len(parts) >= 4:
        return [" ".join(parts[:3]), parts[3]]
    return [str(case_id)]


def draw_panel_frame(draw: ImageDraw.ImageDraw, box: Tuple[int, int, int, int], label: str, title: str):
    x, y, w, h = box
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=PANEL_BG, outline=BORDER, width=BORDER_WIDTH)
    draw_text(draw, (x + 34, y + 38), label, size=62, bold=True, anchor="lm")
    draw_text(draw, (x + w / 2, y + 52), title, size=54, bold=True)


def draw_rmsd_panel(draw: ImageDraw.ImageDraw, box: Tuple[int, int, int, int], metrics: pd.DataFrame):
    draw_panel_frame(draw, box, "A", "Phosphate RMSD")
    x, y, w, h = box
    work = metrics[metrics["variant"].isin(VARIANTS)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    values = [safe_float(v) for v in work["rmsd"]]
    finite_values = [v for v in values if not math.isnan(v)]
    max_y = max(finite_values + [1.0]) * 1.18

    left = x + 235
    top = y + 170
    plot_w = w - 395
    plot_h = h - 510
    group_w = plot_w / max(1, len(cases))
    bar_w = min(50, group_w / (len(VARIANTS) + 1.8))

    draw.line((left, top, left, top + plot_h), fill="#333333", width=5)
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill="#333333", width=5)

    for i in range(6):
        tick = max_y * i / 5
        yy = top + plot_h - (tick / max_y) * plot_h
        draw.line((left - 9, yy, left + plot_w, yy), fill=GRID, width=2)
        draw_text(draw, (left - 22, yy), f"{tick:.1f}", size=32, anchor="rm", fill=MUTED)

    for i, case_id in enumerate(cases):
        cx = left + group_w * i + group_w / 2
        sub = work[work["case_id"] == case_id]
        for j, variant in enumerate(VARIANTS):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = safe_float(row.iloc[0]["rmsd"])
            if math.isnan(value):
                continue
            bx = cx - (len(VARIANTS) * bar_w) / 2 + j * bar_w
            bh = (value / max_y) * plot_h
            by = top + plot_h - bh
            draw.rectangle((bx, by, bx + bar_w - 6, top + plot_h), fill=PANEL_COLORS[variant])
        for k, part in enumerate(wrap_case_label(case_id)):
            draw_text(draw, (cx, top + plot_h + 54 + k * 40), part, size=31, fill="#333333")

    draw_text(draw, (x + 76, top + plot_h / 2), "RMSD (A)", size=36, fill="#333333")

    legend_y = y + h - 88
    legend_x = left + 25
    spacing = 405
    for i, variant in enumerate(VARIANTS):
        lx = legend_x + i * spacing
        draw.rectangle((lx, legend_y - 16, lx + 34, legend_y + 18), fill=PANEL_COLORS[variant])
        draw_text(draw, (lx + 48, legend_y + 1), variant, size=31, anchor="lm", fill="#333333")


def angle_point(cx: float, cy: float, radius: float, angle: float) -> Tuple[float, float]:
    theta = math.radians(float(angle) - 90.0)
    return cx + radius * math.cos(theta), cy + radius * math.sin(theta)


def signed_delta(target: float, reference: float) -> float:
    return ((float(target) - float(reference) + 180.0) % 360.0) - 180.0


def arc_points(cx: float, cy: float, radius: float, reference: float, target: float):
    delta = signed_delta(target, reference)
    steps = max(14, int(abs(delta) / 4.0) + 1)
    return [angle_point(cx, cy, radius, float(reference) + delta * i / steps) for i in range(steps + 1)]


def mix(color: str, other: str = "#ffffff", weight: float = 0.52):
    color = color.lstrip("#")
    other = other.lstrip("#")
    a = tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))
    b = tuple(int(other[i:i + 2], 16) for i in (0, 2, 4))
    return tuple(int(a[i] * (1 - weight) + b[i] * weight) for i in range(3))


def draw_torsion_panel(draw: ImageDraw.ImageDraw, box: Tuple[int, int, int, int], metrics: pd.DataFrame):
    draw_panel_frame(draw, box, "B", "Torsion mode")
    x, y, w, h = box
    work = metrics[metrics["variant"].isin(VARIANTS)].copy()
    cases = list(dict.fromkeys(work["case_id"]))

    cols = 3
    outer_r = 165
    inner_r = 68
    point_outer = 147
    cell_w = (w - 180) / cols
    cell_h = 525
    start_x = x + 90
    start_y = y + 140

    for idx, case_id in enumerate(cases):
        col = idx % cols
        row = idx // cols
        cx = start_x + cell_w * col + cell_w / 2
        cy = start_y + cell_h * row + 168

        for radius in (inner_r, 98, point_outer, outer_r):
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), outline="#e2e2e2", width=3)
        for tick in (-180, -90, 0, 90, 180):
            tx, ty = angle_point(cx, cy, outer_r, tick)
            draw.line((cx, cy, tx, ty), fill="#eeeeee", width=2)
            if tick != -180:
                lx, ly = angle_point(cx, cy, outer_r + 36, tick)
                label = "+/-180" if tick == 180 else str(tick)
                draw_text(draw, (lx, ly), label, size=26, fill=MUTED)
        draw.ellipse((cx - outer_r, cy - outer_r, cx + outer_r, cy + outer_r), outline="#222222", width=5)

        sub = work[work["case_id"] == case_id]
        ref = safe_float(sub.iloc[0]["torsion_ref"]) if not sub.empty else float("nan")
        if not math.isnan(ref):
            rx, ry = angle_point(cx, cy, outer_r + 10, ref)
            draw.line((cx, cy, rx, ry), fill=PANEL_COLORS["Reference"], width=7)
            draw.ellipse((rx - 10, ry - 10, rx + 10, ry + 10), fill=PANEL_COLORS["Reference"])

        for j, variant in enumerate(VARIANTS):
            row_data = sub[sub["variant"] == variant]
            if row_data.empty:
                continue
            value = safe_float(row_data.iloc[0]["torsion_pred"])
            if math.isnan(value):
                continue
            radius = inner_r + j * (point_outer - inner_r) / (len(VARIANTS) - 1)
            color = PANEL_COLORS[variant]
            if not math.isnan(ref):
                draw.line(arc_points(cx, cy, radius, ref, value), fill=mix(color), width=7)
            px, py = angle_point(cx, cy, radius, value)
            draw.ellipse((px - 12, py - 12, px + 12, py + 12), fill=color, outline="#222222", width=2)

        for k, part in enumerate(wrap_case_label(case_id)):
            draw_text(draw, (cx, cy + outer_r + 62 + k * 38), part, size=31, bold=(k == 0))

    legend_items = ["Reference"] + list(VARIANTS)
    legend_y = y + h - 86
    legend_x = x + 260
    spacing = 345
    for i, item in enumerate(legend_items):
        lx = legend_x + i * spacing
        if item == "Reference":
            draw.line((lx, legend_y, lx + 50, legend_y), fill=PANEL_COLORS[item], width=6)
            draw.ellipse((lx + 40, legend_y - 9, lx + 58, legend_y + 9), fill=PANEL_COLORS[item])
        else:
            color = PANEL_COLORS[item]
            draw.line((lx, legend_y, lx + 50, legend_y), fill=mix(color), width=6)
            draw.ellipse((lx + 14, legend_y - 11, lx + 36, legend_y + 11), fill=color, outline="#222222", width=2)
        draw_text(draw, (lx + 66, legend_y), item, size=28, anchor="lm", fill="#333333")


def recovery_fill(recovered: int, total: int) -> str:
    if total <= 0:
        return "#eeeeee"
    ratio = recovered / total
    if ratio >= 0.999:
        return "#cfead1"
    if ratio > 0:
        return "#fff1bd"
    return "#f4c7c3"


def draw_recovery_panel(draw: ImageDraw.ImageDraw, box: Tuple[int, int, int, int], env: pd.DataFrame):
    x, y, w, h = box
    work = env[env["variant"].isin(VARIANTS)].copy()
    cases = list(dict.fromkeys(work["case_id"]))

    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=PANEL_BG, outline=BORDER, width=BORDER_WIDTH)
    draw_text(draw, (x + 34, y + 38), "C", size=62, bold=True, anchor="lm")
    draw_text(draw, (x + w / 2, y + 55), "Basic contact recovery", size=54, bold=True)

    table_pad_x = 70
    left = x + table_pad_x
    right = x + w - table_pad_x
    table_w = right - left
    top = y + 160
    row_h = 91
    col_case = 600
    col_crystal = 190
    tool_gap = 38
    variant_area_w = table_w - col_case - col_crystal
    cell_w = (variant_area_w - (len(VARIANTS) - 1) * tool_gap) / max(1, len(VARIANTS))
    variant_x0 = left + col_case + col_crystal

    draw.line((left, top + 36, right, top + 36), fill="#dddddd", width=3)
    draw_text(draw, (left, top), "Case", size=36, anchor="lm", bold=True)
    draw_text(draw, (left + col_case + 42, top), "Crystal", size=36, anchor="lm", bold=True)
    for i, variant in enumerate(VARIANTS):
        hx = variant_x0 + i * (cell_w + tool_gap) + cell_w / 2
        draw_text(draw, (hx, top), variant, size=34, bold=True)

    for i, case_id in enumerate(cases):
        yy = top + 84 + i * row_h
        sub = work[work["case_id"] == case_id]
        first = sub.iloc[0]
        total = int(first.get("phospho_crystal_basic_contact_count_4A", 0))
        if i % 2:
            draw.rectangle((left - 30, yy - 42, right + 30, yy + 42), fill="#fafafa")
        draw_text(draw, (left, yy), case_id, size=34, anchor="lm", fill="#333333")
        draw_text(draw, (left + col_case + 82, yy), str(total), size=36, fill="#333333")
        for j, variant in enumerate(VARIANTS):
            row = sub[sub["variant"] == variant]
            cell_x = variant_x0 + j * (cell_w + tool_gap)
            label = "NA"
            fill = "#eeeeee"
            if not row.empty:
                recovered = int(row.iloc[0].get("recovered_crystal_basic_contact_count_4A", 0))
                pct = int(round(100 * recovered / total)) if total else 0
                label = f"{recovered}/{total} ({pct}%)"
                fill = recovery_fill(recovered, total)
            draw.rounded_rectangle((cell_x, yy - 34, cell_x + cell_w, yy + 34), radius=8, fill=fill, outline="#dddddd", width=2)
            draw_text(draw, (cell_x + cell_w / 2, yy), label, size=34)


def load_rgb(path: Path) -> Image.Image:
    image = Image.open(path)
    if image.mode == "RGBA":
        bg = Image.new("RGBA", image.size, "white")
        bg.alpha_composite(image)
        return bg.convert("RGB")
    return image.convert("RGB")


def paste_pymol_panel(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    image_path: Path,
    label: str,
    title: str,
):
    x, y, w, h = box
    image = load_rgb(image_path)
    fitted = ImageOps.fit(image, (w, h), method=Image.Resampling.LANCZOS, centering=(0.50, 0.48))
    canvas.paste(fitted, (x, y))
    draw_text(draw, (x + 34, y + 38), label, size=62, bold=True, anchor="lm")

    title_font_size = 44
    title_w = text_width(draw, title, title_font_size, bold=True)
    draw.rectangle((x, y + h - 74, x + title_w + 52, y + h), fill="white")
    draw_text(draw, (x + 26, y + h - 44), title, size=title_font_size, bold=True, anchor="lm")
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, outline=BORDER, width=BORDER_WIDTH)


def assert_inputs_exist(paths: Sequence[Path]):
    missing = [path for path in paths if not path.exists()]
    if missing:
        joined = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing required input files:\n{joined}")


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


def main() -> None:
    metrics_path = OUT_DIR / "case_study_metrics_long.tsv"
    env_path = OUT_DIR / "case_study_crystal_basic_contact_recovery.tsv"
    pymol_panels = [
        (ROOT / "CDK2_PF1.png", "D", "PhosphoFill top1"),
        (ROOT / "CDK2_Pytm-def.png", "E", "PyTMs default"),
        (ROOT / "CDK2_Pytm-opt.png", "F", "PyTMs optimized"),
        (ROOT / "CDK2_PTM-psi.png", "G", "PTM-Psi"),
    ]
    assert_inputs_exist([metrics_path, env_path] + [path for path, _, _ in pymol_panels])

    metrics = pd.read_csv(metrics_path, sep="\t")
    env = pd.read_csv(env_path, sep="\t")

    image = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)

    margin_x = 75
    top_y = 58
    top_h = 1418
    top_gap = 70
    top_w = (WIDTH - 2 * margin_x - top_gap) // 2
    draw_rmsd_panel(draw, (margin_x, top_y, top_w, top_h), metrics)
    draw_torsion_panel(draw, (margin_x + top_w + top_gap, top_y, top_w, top_h), metrics)

    recovery_y = 1545
    recovery_h = 865
    draw_recovery_panel(draw, (margin_x, recovery_y, WIDTH - 2 * margin_x, recovery_h), env)

    draw_text(
        draw,
        (WIDTH / 2, 2488),
        "Phosphate contact recovery compared between predictions",
        size=54,
        bold=True,
    )

    pymol_y = 2572
    pymol_h = 965
    pymol_x = 55
    gap = 30
    pymol_w = (WIDTH - 2 * pymol_x - 3 * gap) // 4
    for i, (path, label, title) in enumerate(pymol_panels):
        px = pymol_x + i * (pymol_w + gap)
        paste_pymol_panel(image, draw, (px, pymol_y, pymol_w, pymol_h), path, label, title)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    image.save(OUTPUT_PNG, dpi=(DPI, DPI))
    image.save(OUTPUT_TIF, dpi=(DPI, DPI), compression="tiff_lzw")
    save_pdf(image, OUTPUT_PDF)
    print(OUTPUT_PNG)
    print(OUTPUT_TIF)
    print(OUTPUT_PDF)
    print(image.size)


if __name__ == "__main__":
    main()
