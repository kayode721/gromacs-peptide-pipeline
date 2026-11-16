#!/usr/bin/env python3
"""
Plot a GROMACS DSSP .dat file as a publication-style secondary-structure map.

Examples:
  python input/plot_dssp_states.py saa1-27nh2.dat --duration-ns 100
  python input/plot_dssp_states.py saa1-27nh2.dat --dt-ns 0.002 --subsample 10
"""

import argparse
import csv
import html
from pathlib import Path


DETAILED_STATES = [
    ("~", "Coil", "#ffffff"),
    ("E", "Beta strand", "#d73027"),
    ("B", "Beta bridge", "#000000"),
    ("S", "Bend", "#1a9850"),
    ("T", "Turn", "#fee08b"),
    ("H", "Alpha helix", "#7b3294"),
    ("G", "3-10 helix", "#4575b4"),
    ("I", "Pi helix", "#74add1"),
    ("P", "PPII", "#66c2a5"),
    ("=", "Break", "#bdbdbd"),
]

SIMPLE_STATES = [
    ("coil", "Coil/other", "#ffffff"),
    ("beta", "Beta strand/bridge", "#d73027"),
    ("turn", "Turn/bend", "#fee08b"),
    ("helix", "Helix", "#7b3294"),
]


def parse_residue_selection(selection):
    if not selection:
        return None

    residues = set()
    for part in selection.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            residues.update(range(int(start), int(end) + 1))
        else:
            residues.add(int(part))

    return sorted(residues)


def read_dssp_dat(path):
    lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = "".join(raw.split())
            if line:
                lines.append(line)

    if not lines:
        raise ValueError(f"No DSSP rows found in {path}")

    width = len(lines[0])
    bad = [i + 1 for i, line in enumerate(lines) if len(line) != width]
    if bad:
        raise ValueError(f"DSSP rows are not all the same length. First bad row: {bad[0]}")

    return lines


def state_definitions(mode):
    if mode == "simple":
        return SIMPLE_STATES
    return DETAILED_STATES


def detailed_code(letter):
    letter = letter.upper()
    if letter in {"H", "G", "I", "E", "B", "P", "T", "S", "="}:
        return letter
    return "~"


def simple_code(letter):
    letter = letter.upper()
    if letter in {"H", "G", "I"}:
        return "helix"
    if letter in {"E", "B"}:
        return "beta"
    if letter in {"T", "S"}:
        return "turn"
    return "coil"


def encode_rows(lines, mode):
    states = state_definitions(mode)
    lookup = {code: i for i, (code, _label, _color) in enumerate(states)}
    coder = simple_code if mode == "simple" else detailed_code
    return [[lookup[coder(letter)] for letter in line] for line in lines]


def transpose(rows):
    return [list(col) for col in zip(*rows)]


def make_time_axis(raw_frame_count, kept_indices, dt_ns, duration_ns):
    if duration_ns is not None:
        if raw_frame_count <= 1:
            dt_ns = 0.0
        else:
            dt_ns = duration_ns / (raw_frame_count - 1)
    elif dt_ns is None:
        return [float(i) for i in kept_indices], "Frame index"

    return [i * dt_ns for i in kept_indices], "Time (ns)"


def write_summary_csv(path, lines, residue_numbers, mode):
    states = state_definitions(mode)
    counts = [{code: 0 for code, _label, _color in states} for _ in residue_numbers]
    coder = simple_code if mode == "simple" else detailed_code

    for line in lines:
        for pos, letter in enumerate(line):
            counts[pos][coder(letter)] += 1

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["residue", *[label for _code, label, _color in states]])
        for residue, residue_counts in zip(residue_numbers, counts):
            total = sum(residue_counts.values()) or 1
            writer.writerow([
                residue,
                *[round(100.0 * residue_counts[code] / total, 3) for code, _label, _color in states],
            ])


def fmt_number(value):
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def write_svg(data, xvals, x_label, residue_numbers, states, title, out):
    width = 1200
    row_h = 12
    left = 76
    right = 24
    top = 44 if title else 24
    bottom = 108
    plot_w = width - left - right
    plot_h = max(1, len(residue_numbers)) * row_h
    height = top + plot_h + bottom
    nframes = len(xvals)
    colors = [color for _code, _label, color in states]

    def x_for_index(index):
        return left + (index / max(1, nframes)) * plot_w

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>',
        'text { font-family: Arial, Helvetica, sans-serif; fill: #222; }',
        '.axis { stroke: #222; stroke-width: 1; }',
        '.tick { stroke: #222; stroke-width: 1; }',
        '.grid { stroke: #e5e5e5; stroke-width: 1; }',
        '</style>',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
    ]

    if title:
        lines.append(
            f'<text x="{width / 2:.1f}" y="24" text-anchor="middle" font-size="18" font-weight="700">'
            f'{html.escape(title)}</text>'
        )

    lines.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f7f7f7"/>')

    for row_i, row in enumerate(data):
        y = top + (len(data) - 1 - row_i) * row_h
        start = 0
        current = row[0]
        for idx in range(1, len(row) + 1):
            if idx == len(row) or row[idx] != current:
                x = x_for_index(start)
                rect_w = max(0.4, x_for_index(idx) - x)
                lines.append(
                    f'<rect x="{x:.3f}" y="{y:.3f}" width="{rect_w:.3f}" height="{row_h}" '
                    f'fill="{colors[current]}"/>'
                )
                if idx < len(row):
                    start = idx
                    current = row[idx]

    # Axes and ticks.
    axis_y = top + plot_h
    lines.append(f'<line class="axis" x1="{left}" y1="{axis_y}" x2="{left + plot_w}" y2="{axis_y}"/>')
    lines.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{axis_y}"/>')

    tick_count = 10
    for t in range(tick_count + 1):
        idx = round(t * (nframes - 1) / tick_count) if nframes > 1 else 0
        x = x_for_index(idx)
        label = fmt_number(xvals[idx])
        lines.append(f'<line class="grid" x1="{x:.3f}" y1="{top}" x2="{x:.3f}" y2="{axis_y}"/>')
        lines.append(f'<line class="tick" x1="{x:.3f}" y1="{axis_y}" x2="{x:.3f}" y2="{axis_y + 5}"/>')
        lines.append(f'<text x="{x:.3f}" y="{axis_y + 22}" text-anchor="middle" font-size="12">{label}</text>')

    if len(residue_numbers) <= 40:
        y_ticks = list(range(len(residue_numbers)))
    else:
        step = max(1, len(residue_numbers) // 20)
        y_ticks = list(range(0, len(residue_numbers), step))

    for row_i in y_ticks:
        residue = residue_numbers[row_i]
        y = top + (len(data) - 1 - row_i) * row_h + row_h / 2
        lines.append(f'<line class="tick" x1="{left - 5}" y1="{y:.3f}" x2="{left}" y2="{y:.3f}"/>')
        lines.append(
            f'<text x="{left - 9}" y="{y + 4:.3f}" text-anchor="end" font-size="11">{residue}</text>'
        )

    lines.append(
        f'<text x="{left + plot_w / 2:.1f}" y="{height - 44}" text-anchor="middle" font-size="14">'
        f'{html.escape(x_label)}</text>'
    )
    lines.append(
        f'<text x="18" y="{top + plot_h / 2:.1f}" text-anchor="middle" font-size="14" '
        f'transform="rotate(-90 18 {top + plot_h / 2:.1f})">Residue number</text>'
    )

    legend_y = height - 26
    legend_x = left
    for _code, label, color in states:
        escaped = html.escape(label)
        lines.append(
            f'<rect x="{legend_x}" y="{legend_y - 12}" width="12" height="12" '
            f'fill="{color}" stroke="black" stroke-width="0.8"/>'
        )
        lines.append(f'<text x="{legend_x + 18}" y="{legend_y - 2}" font-size="12">{escaped}</text>')
        legend_x += 18 + max(68, 7 * len(label))

    lines.append("</svg>")

    with open(out, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print(f"Saved figure to {out}")


def plot_dssp(data, xvals, x_label, residue_numbers, states, title, out, dpi, show):
    if str(out).lower().endswith(".svg"):
        write_svg(data, xvals, x_label, residue_numbers, states, title, out)
        return

    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm, ListedColormap
        from matplotlib.patches import Rectangle
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This script needs matplotlib. Run it with the same Python environment "
            "you used for input/plot_dssp.py, install matplotlib, or use an .svg output file."
        ) from exc

    colors = [color for _code, _label, color in states]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm([i - 0.5 for i in range(len(states) + 1)], cmap.N)

    height = max(3.0, min(7.0, 0.22 * len(residue_numbers) + 1.5))
    fig, ax = plt.subplots(figsize=(12, height))
    ax.imshow(
        data,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        origin="lower",
        extent=[
            xvals[0],
            xvals[-1] if len(xvals) > 1 else xvals[0] + 1,
            residue_numbers[0] - 0.5,
            residue_numbers[-1] + 0.5,
        ],
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel("Residue number")
    if title:
        ax.set_title(title)

    if len(xvals) > 1:
        step_x = max(1, len(xvals) // 10)
        ax.set_xticks([xvals[i] for i in range(0, len(xvals), step_x)])
    if len(residue_numbers) <= 40:
        ax.set_yticks(residue_numbers)
    else:
        step_y = max(1, len(residue_numbers) // 20)
        ax.set_yticks(residue_numbers[::step_y])

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="black", label=label)
        for _code, label, color in states
    ]
    ax.legend(
        handles=legend_handles,
        bbox_to_anchor=(0.5, -0.22),
        loc="upper center",
        ncol=min(len(states), 5),
        frameon=False,
        handleheight=1.2,
        handlelength=1.2,
    )

    fig.tight_layout(rect=[0, 0.12, 1, 1])
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"Saved figure to {out}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Make a labeled DSSP secondary-structure map with time on the x-axis."
    )
    parser.add_argument("dat_file", help="GROMACS DSSP .dat file")
    parser.add_argument("--dt-ns", type=float, help="Time spacing between DSSP frames in ns")
    parser.add_argument("--duration-ns", type=float, help="Total trajectory duration in ns")
    parser.add_argument("--res-offset", type=int, default=1, help="First residue number")
    parser.add_argument("--residues", help="Residue range, e.g. 1-10,15,20-27")
    parser.add_argument("--subsample", type=int, default=1, help="Keep every Nth DSSP frame")
    parser.add_argument(
        "--mode",
        choices=["detailed", "simple"],
        default="detailed",
        help="detailed shows DSSP states; simple groups helix/beta/turn/coil",
    )
    parser.add_argument("--title", default="", help="Optional plot title")
    parser.add_argument("--out", help="Output image path, usually .svg or .png")
    parser.add_argument("--summary-csv", help="Output per-residue percent occupancy CSV")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show", action="store_true", help="Open an interactive plot window")
    args = parser.parse_args()

    if args.subsample < 1:
        raise ValueError("--subsample must be >= 1")
    if args.dt_ns is not None and args.duration_ns is not None:
        raise ValueError("Use either --dt-ns or --duration-ns, not both")

    dat_path = Path(args.dat_file)
    raw_lines = read_dssp_dat(dat_path)
    raw_frame_count = len(raw_lines)
    kept_indices = list(range(0, raw_frame_count, args.subsample))
    lines = [raw_lines[i] for i in kept_indices]

    residue_numbers = list(range(args.res_offset, args.res_offset + len(lines[0])))
    selection = parse_residue_selection(args.residues)
    if selection is not None:
        positions = [i for i, residue in enumerate(residue_numbers) if residue in selection]
        if not positions:
            raise ValueError("None of the selected residues are present")
        residue_numbers = [residue_numbers[i] for i in positions]
        lines = ["".join(line[i] for i in positions) for line in lines]

    encoded = encode_rows(lines, args.mode)
    data = transpose(encoded)
    xvals, x_label = make_time_axis(raw_frame_count, kept_indices, args.dt_ns, args.duration_ns)
    states = state_definitions(args.mode)

    out = args.out or str(dat_path.with_name(f"{dat_path.stem}_dssp_states.svg"))
    plot_dssp(data, xvals, x_label, residue_numbers, states, args.title, out, args.dpi, args.show)

    summary_csv = args.summary_csv
    if summary_csv is None:
        summary_csv = str(dat_path.with_name(f"{dat_path.stem}_dssp_summary.csv"))
    write_summary_csv(summary_csv, lines, residue_numbers, args.mode)
    print(f"Saved summary to {summary_csv}")


if __name__ == "__main__":
    main()
