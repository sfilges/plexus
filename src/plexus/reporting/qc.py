"""Panel QC report generation (REPT-01)."""

from __future__ import annotations

import csv
import re
import statistics
from itertools import combinations
from typing import TYPE_CHECKING

from plexus.aligner import PrimerDimerPredictor

if TYPE_CHECKING:
    from plexus.designer.multiplexpanel import Junction


def _tailed(tail: str, seq: str) -> str:
    """Prepend tail to sequence, replacing N with A for NN-model compatibility."""
    return tail.replace("N", "A") + seq


def generate_panel_qc(
    junctions: list[Junction],
    *,
    gc_high_threshold: float = 70.0,
    gc_low_threshold: float = 30.0,
    homopolymer_min_run: int = 4,
    dimer_threshold: float = 0.0,
    forward_tail: str = "",
    reverse_tail: str = "",
) -> dict:
    """Generate panel QC metrics for the selected primer pairs."""
    # Build working lists
    selected = []  # list of (junction_name, PrimerPair)
    for j in junctions:
        for pair in j.primer_pairs:
            if pair.selected:
                selected.append((j.name, pair))
                break

    all_primers = []  # (junction_name, direction, Primer)
    for jname, pair in selected:
        all_primers.append((jname, "forward", pair.forward))
        all_primers.append((jname, "reverse", pair.reverse))

    # Tm distribution
    tms = [p.tm for _, _, p in all_primers]
    tm_distribution = {
        "mean": round(statistics.mean(tms), 2) if tms else None,
        "std": round(statistics.stdev(tms), 2) if len(tms) >= 2 else None,
        "min": round(min(tms), 2) if tms else None,
        "max": round(max(tms), 2) if tms else None,
        "per_primer": [
            {"junction": jname, "direction": d, "name": p.name, "tm": round(p.tm, 2)}
            for jname, d, p in all_primers
        ],
    }

    # Sequence flags
    _hp_re = re.compile(
        r"([ACGT])\1{" + str(homopolymer_min_run - 1) + r",}", re.IGNORECASE
    )
    flagged_primers = []
    high_gc_count = low_gc_count = homopolymer_count = 0
    for jname, direction, p in all_primers:
        flags = []
        if p.gc > gc_high_threshold:
            flags.append("high_gc")
            high_gc_count += 1
        if p.gc < gc_low_threshold:
            flags.append("low_gc")
            low_gc_count += 1
        if _hp_re.search(p.seq):
            flags.append("homopolymer")
            homopolymer_count += 1
        if flags:
            flagged_primers.append(
                {
                    "junction": jname,
                    "direction": direction,
                    "name": p.name,
                    "gc": round(p.gc, 1),
                    "sequence": p.seq,
                    "flags": flags,
                }
            )
    sequence_flags = {
        "gc_high_threshold": gc_high_threshold,
        "gc_low_threshold": gc_low_threshold,
        "homopolymer_min_run": homopolymer_min_run,
        "high_gc_count": high_gc_count,
        "low_gc_count": low_gc_count,
        "homopolymer_count": homopolymer_count,
        "flagged_primers": flagged_primers,
    }

    # --- Primer-level dimer matrix ---
    # Build list of (label, tailed_sequence) for every primer in the panel.
    primer_entries = []  # [(label, tailed_seq), ...]
    primer_labels = []
    for jname, pair in selected:
        fwd_label = f"{jname}_forward"
        rev_label = f"{jname}_reverse"
        primer_labels.append(fwd_label)
        primer_labels.append(rev_label)
        primer_entries.append((fwd_label, _tailed(forward_tail, pair.forward.seq)))
        primer_entries.append((rev_label, _tailed(reverse_tail, pair.reverse.seq)))

    predictor = PrimerDimerPredictor()
    primer_matrix: dict[str, dict[str, float]] = {}

    for (label_a, seq_a), (label_b, seq_b) in combinations(primer_entries, 2):
        predictor.set_primers(seq_a, seq_b, label_a, label_b)
        predictor.align()
        score = round(predictor.score or 0.0, 4)
        primer_matrix.setdefault(label_a, {})[label_b] = score
        primer_matrix.setdefault(label_b, {})[label_a] = score

    primer_dimer_matrix = {
        "primer_labels": primer_labels,
        "dimer_threshold": dimer_threshold,
        "matrix": primer_matrix,
    }

    # --- Junction-level cross-reactivity matrix (derived) ---
    junction_matrix: dict[str, dict] = {}
    for (jname_a, _pair_a), (jname_b, _pair_b) in combinations(selected, 2):
        scores = []
        for d_a in ("forward", "reverse"):
            for d_b in ("forward", "reverse"):
                la = f"{jname_a}_{d_a}"
                lb = f"{jname_b}_{d_b}"
                scores.append(primer_matrix.get(la, {}).get(lb, 0.0))
        cell = {
            "min_dimer_score": round(min(scores), 4),
            "interaction_count": sum(1 for s in scores if s < dimer_threshold),
        }
        junction_matrix.setdefault(jname_a, {})[jname_b] = cell
        junction_matrix.setdefault(jname_b, {})[jname_a] = cell

    cross_reactivity_matrix = {
        "dimer_threshold": dimer_threshold,
        "matrix": junction_matrix,
    }

    return {
        "tm_distribution": tm_distribution,
        "sequence_flags": sequence_flags,
        "cross_reactivity_matrix": cross_reactivity_matrix,
        "primer_dimer_matrix": primer_dimer_matrix,
    }


def save_primer_dimer_matrix_csv(
    primer_dimer_data: dict,
    file_path: str,
) -> None:
    """Write the primer-level dimer matrix as a symmetric CSV."""
    labels = primer_dimer_data["primer_labels"]
    matrix = primer_dimer_data["matrix"]

    with open(file_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + labels)
        for label_a in labels:
            row = [label_a]
            for label_b in labels:
                if label_a == label_b:
                    row.append("")
                else:
                    row.append(matrix.get(label_a, {}).get(label_b, ""))
            writer.writerow(row)
