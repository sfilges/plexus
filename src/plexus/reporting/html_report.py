"""HTML QC report generation (REPT-02).

Transforms panel_qc.json and related pipeline outputs into a self-contained,
interactive HTML report with Plotly.js charts.
"""

from __future__ import annotations

import csv
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from plexus.version import __version__

_TEMPLATES_DIR = Path(__file__).parent / "templates"
_PLOTLY_JS_GZ_PATH = _TEMPLATES_DIR / "plotly.min.js.gz"


def _read_json(path: Path) -> dict | None:
    if path.is_file():
        with path.open() as f:
            return json.load(f)
    return None


def _read_csv(path: Path) -> list[dict] | None:
    if path.is_file():
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            rows = []
            for row in reader:
                # Convert numeric fields
                for key in row:
                    val = row[key]
                    if val == "":
                        continue
                    try:
                        row[key] = int(val)
                    except (ValueError, TypeError):
                        try:
                            row[key] = float(val)
                        except (ValueError, TypeError):
                            pass
                rows.append(row)
            return rows if rows else None
    return None


def _load_plotly_js() -> str:
    with gzip.open(_PLOTLY_JS_GZ_PATH, "rt", encoding="utf-8") as f:
        return f.read()


def _render_report(
    qc_data: dict,
    *,
    panel_name: str = "Panel",
    summary_data: dict | None = None,
    selected_pairs: list[dict] | None = None,
    off_targets: list[dict] | None = None,
    failed_junctions: list[dict] | None = None,
    provenance_data: dict | None = None,
    top_panels: list[dict] | None = None,
) -> str:
    """Render the HTML report string from data."""
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape=False,
    )
    template = env.get_template("panel_report.html.j2")

    num_failed = len(failed_junctions) if failed_junctions else 0

    return template.render(
        panel_name=panel_name,
        report_date=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        qc=qc_data,
        qc_json=json.dumps(qc_data),
        summary=summary_data,
        num_failed=num_failed,
        provenance=provenance_data,
        selected_pairs_data=selected_pairs,
        selected_pairs_json=json.dumps(selected_pairs) if selected_pairs else "[]",
        off_targets_data=off_targets,
        failed_junctions_data=failed_junctions,
        top_panels_data=top_panels,
        plexus_version=__version__,
        plotly_js=_load_plotly_js(),
    )


def generate_html_report(
    output_dir: Path,
    *,
    panel_name: str = "Panel",
) -> Path:
    """Generate an HTML QC report from pipeline output files.

    Reads from output_dir: panel_qc.json (required), panel_summary.json,
    selected_multiplex.csv, off_targets.csv, failed_junctions.csv,
    provenance.json, top_panels.csv — all optional except panel_qc.json.

    Returns the path to the generated panel_report.html.
    """
    output_dir = Path(output_dir)

    qc_data = _read_json(output_dir / "panel_qc.json")
    if qc_data is None:
        raise FileNotFoundError(f"panel_qc.json not found in {output_dir}")

    summary_data = _read_json(output_dir / "panel_summary.json")
    provenance_data = _read_json(output_dir / "provenance.json")
    selected_pairs = _read_csv(output_dir / "selected_multiplex.csv")
    off_targets = _read_csv(output_dir / "off_targets.csv")
    failed_junctions = _read_csv(output_dir / "failed_junctions.csv")
    top_panels = _read_csv(output_dir / "top_panels.csv")

    # Use panel name from summary if available
    if summary_data and summary_data.get("panel_name"):
        panel_name = summary_data["panel_name"]

    html = _render_report(
        qc_data,
        panel_name=panel_name,
        summary_data=summary_data,
        selected_pairs=selected_pairs,
        off_targets=off_targets,
        failed_junctions=failed_junctions,
        provenance_data=provenance_data,
        top_panels=top_panels,
    )

    html_path = output_dir / "panel_report.html"
    html_path.write_text(html, encoding="utf-8")
    return html_path


def generate_html_report_from_data(
    qc_data: dict,
    summary_data: dict | None,
    selected_pairs_csv: str | None,
    off_targets_csv: str | None,
    failed_junctions_csv: str | None,
    output_path: Path,
    *,
    panel_name: str = "Panel",
    provenance_data: dict | None = None,
    top_panels_csv: str | None = None,
) -> Path:
    """Generate HTML report from in-memory data (for pipeline integration).

    CSV arguments are raw CSV strings; they are parsed internally.
    """
    import io

    def _parse_csv_string(csv_str: str | None) -> list[dict] | None:
        if not csv_str:
            return None
        reader = csv.DictReader(io.StringIO(csv_str))
        rows = []
        for row in reader:
            for key in row:
                val = row[key]
                if val == "":
                    continue
                try:
                    row[key] = int(val)
                except (ValueError, TypeError):
                    try:
                        row[key] = float(val)
                    except (ValueError, TypeError):
                        pass
            rows.append(row)
        return rows if rows else None

    selected_pairs = _parse_csv_string(selected_pairs_csv)
    off_targets = _parse_csv_string(off_targets_csv)
    failed_junctions = _parse_csv_string(failed_junctions_csv)
    top_panels = _parse_csv_string(top_panels_csv)

    if summary_data and summary_data.get("panel_name"):
        panel_name = summary_data["panel_name"]

    html = _render_report(
        qc_data,
        panel_name=panel_name,
        summary_data=summary_data,
        selected_pairs=selected_pairs,
        off_targets=off_targets,
        failed_junctions=failed_junctions,
        provenance_data=provenance_data,
        top_panels=top_panels,
    )

    output_path = Path(output_path)
    output_path.write_text(html, encoding="utf-8")
    return output_path
