#!/usr/bin/env python3
"""Create small test fixtures from real reference data for integration tests.

Extracts tiny regions around two junctions and writes:
  tests/data/fixtures.fa       - 2 contigs x 1201bp
  tests/data/fixtures.fa.fai   - FASTA index
  tests/data/junctions.csv     - 3 junctions (2 merge on chr22, 1 on chr7)
  tests/data/snps.vcf.gz       - gnomAD variants remapped to small contigs
  tests/data/snps.vcf.gz.tbi   - tabix index

Requires real data files in data/:
  data/hg38/hg38.fa            - reference genome
  data/af-only-gnomad.hg38.vcf.gz  - gnomAD VCF
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pysam

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "tests" / "data"

# Real genomic coordinates for two junctions
REGIONS = [
    {
        "name": "CLTCL1_p.R354H",
        "chrom": "chr22",
        "position": 19234615,  # 1-based
    },
    {
        "name": "BRAF_p.L485W",
        "chrom": "chr7",
        "position": 140778054,  # 1-based
    },
]

CONTIG_HALF = 600  # bases each side of junction -> 1201bp contig
FIXTURE_POSITION = 601  # 1-based position of junction in small contig


def main() -> None:
    fasta_path = PROJECT_ROOT / "data" / "hg38" / "hg38.fa"
    vcf_path = PROJECT_ROOT / "data" / "af-only-gnomad.hg38.vcf.gz"

    for p in (fasta_path, vcf_path):
        if not p.exists():
            print(f"ERROR: {p} not found", file=sys.stderr)
            sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Extract FASTA contigs
    # ------------------------------------------------------------------
    fa_out = OUTPUT_DIR / "fixtures.fa"
    print(f"Extracting FASTA contigs to {fa_out} ...")

    contig_info: dict[str, dict] = {}  # chrom -> {orig_start, orig_end, seq}

    with pysam.FastaFile(str(fasta_path)) as fasta:
        with open(fa_out, "w") as fh:
            for region in REGIONS:
                chrom = region["chrom"]
                pos = region["position"]

                # 1-based inclusive range in original genome
                orig_start = pos - CONTIG_HALF  # 1-based start
                orig_end = pos + CONTIG_HALF  # 1-based end

                # pysam uses 0-based half-open
                seq = fasta.fetch(chrom, orig_start - 1, orig_end)

                assert (
                    len(seq) == 2 * CONTIG_HALF + 1
                ), f"Expected {2 * CONTIG_HALF + 1}bp, got {len(seq)} for {chrom}"

                fh.write(f">{chrom}\n")
                # Write in 80-char lines
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")

                contig_info[chrom] = {
                    "orig_start": orig_start,
                    "orig_end": orig_end,
                    "seq_len": len(seq),
                }

    # Index the FASTA
    pysam.faidx(str(fa_out))
    print(f"  Indexed: {fa_out}.fai")

    for chrom, info in contig_info.items():
        print(
            f"  {chrom}: {info['seq_len']}bp "
            f"(orig {info['orig_start']}-{info['orig_end']})"
        )

    # ------------------------------------------------------------------
    # 2. Extract and remap VCF variants
    # ------------------------------------------------------------------
    vcf_out = OUTPUT_DIR / "snps.vcf.gz"
    print(f"\nExtracting VCF variants to {vcf_out} ...")

    total_variants = 0
    per_chrom_variants: dict[str, int] = {}

    with pysam.VariantFile(str(vcf_path)) as vcf_in:
        # Build new header with small contigs
        new_header = pysam.VariantHeader()
        for sample in vcf_in.header.samples:
            new_header.add_sample(sample)

        # Copy INFO fields from original header
        for record in vcf_in.header.records:
            if record.type == "INFO":
                new_header.add_record(record)
            elif record.type == "FILTER":
                new_header.add_record(record)
            elif record.type == "FORMAT":
                new_header.add_record(record)

        # Add contig lines for our small contigs
        for chrom, info in contig_info.items():
            new_header.add_line(f"##contig=<ID={chrom},length={info['seq_len']}>")

        with pysam.VariantFile(str(vcf_out), "wz", header=new_header) as vcf_out_fh:
            for region in REGIONS:
                chrom = region["chrom"]
                info = contig_info[chrom]
                offset = info["orig_start"] - 1  # 0-based offset

                # Fetch variants in the original region (0-based)
                count = 0
                for record in vcf_in.fetch(
                    chrom, info["orig_start"] - 1, info["orig_end"]
                ):
                    new_pos = record.pos - offset  # remap to small contig

                    if new_pos < 1 or new_pos > info["seq_len"]:
                        continue

                    new_rec = vcf_out_fh.new_record()
                    new_rec.contig = chrom
                    new_rec.pos = new_pos
                    new_rec.id = record.id
                    new_rec.alleles = record.alleles
                    new_rec.qual = record.qual

                    # Copy INFO fields
                    for key in record.info:
                        try:
                            new_rec.info[key] = record.info[key]
                        except (KeyError, TypeError):
                            pass

                    # Copy filter
                    if record.filter:
                        for f in record.filter:
                            if f != "PASS":
                                try:
                                    new_rec.filter.add(f)
                                except KeyError:
                                    pass

                    vcf_out_fh.write(new_rec)
                    count += 1

                per_chrom_variants[chrom] = count
                total_variants += count

    # Tabix index
    pysam.tabix_index(str(vcf_out), preset="vcf", force=True)
    print(f"  Indexed: {vcf_out}.tbi")
    print(f"  Total variants: {total_variants}")
    for chrom, count in per_chrom_variants.items():
        print(f"    {chrom}: {count} variants")

    # ------------------------------------------------------------------
    # 3. Write junctions CSV
    # ------------------------------------------------------------------
    csv_out = OUTPUT_DIR / "junctions.csv"
    print(f"\nWriting junctions CSV to {csv_out} ...")

    rows = [
        # Two CLTCL1 entries on chr22 (will merge)
        {
            "Name": "CLTCL1_p.R354H",
            "Chrom": "chr22",
            "Five_Prime_Coordinate": FIXTURE_POSITION,
            "Three_Prime_Coordinate": FIXTURE_POSITION,
        },
        {
            "Name": "CLTCL1_p.R354H",
            "Chrom": "chr22",
            "Five_Prime_Coordinate": FIXTURE_POSITION,
            "Three_Prime_Coordinate": FIXTURE_POSITION,
        },
        # One BRAF entry on chr7 (separate chromosome)
        {
            "Name": "BRAF_p.L485W",
            "Chrom": "chr7",
            "Five_Prime_Coordinate": FIXTURE_POSITION,
            "Three_Prime_Coordinate": FIXTURE_POSITION,
        },
    ]

    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "Name",
                "Chrom",
                "Five_Prime_Coordinate",
                "Three_Prime_Coordinate",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"  {len(rows)} junctions written")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n=== Fixture Summary ===")
    for f in sorted(OUTPUT_DIR.iterdir()):
        size = f.stat().st_size
        print(f"  {f.name:25s}  {size:>8,} bytes")
    print("Done!")


if __name__ == "__main__":
    main()
