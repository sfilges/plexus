# Design Plan: Clinical-Grade In-Silico PCR (isPCR)

## 1. Overview

The current Plexus specificity check identifies where individual primers bind (via BLAST) but does not explicitly model the formation of **amplicons**. In a clinical context (e.g., ctDNA panels), it is critical to identify not just non-specific binding, but non-specific *amplification*. This document outlines a plan to implement a dedicated in-silico PCR (isPCR) module.

## 2. Core Requirements

- **Pair Identification**: Identify pairs of primers (Forward/Reverse) that bind on the same contig in the correct orientation within a specific distance (e.g., < 5000 bp).
- **Cross-Amplification**: Detect amplicons formed by primers from different targets (e.g., Target A Forward + Target B Reverse).
- **Thermodynamic Validation**: Use physical models to predict if a binding site will actually lead to polymerase extension.

## 3. Methodology Comparison

### 3.1 Thermodynamic Model (Johnston/Plexus)

- **Mechanism**: Ungapped alignment, heavy 3'-end bonus.
- **Pros**: Specifically tuned for PCR "extensibility." Excellent at penalizing 3' mismatches that stop polymerization.
- **Cons**: Heuristic; ignores gaps and internal stability (mismatches in the middle of the primer).

### 3.2 $\Delta G$ Model (`ntthal` / Primer3)

- **Mechanism**: Full nearest-neighbor thermodynamic model (SantaLucia 1998).
- **Pros**: Physically accurate; handles gapped alignments and salt/buffer conditions.
- **Cons**: May over-penalize strong internal binding that cannot extend due to a 3' mismatch.

### 3.3 The Hybrid Strategy (Recommended)

We will use a multi-stage filtering approach:

1. **Discovery (BLAST)**: Find all genomic locations with significant homology.
2. **Binding Strength (`ntthal`)**: Calculate the actual $\Delta G$ of binding for every BLAST hit.
3. **Genomic Pairing**: Identify $F/R$ hits within 5kb.
4. **Extension Risk (Johnston)**: Apply the 3'-end stability bonus to the identified pairs to classify the functional risk of amplification.

## 4. Proposed Implementation Path

### Phase 1: Enhance `BlastRunner`

- Update BLAST parameters to be more sensitive (lower `word_size`, higher `e-value`) to ensure weak binding sites that could pair up are not missed.

### Phase 2: Complete `AmpliconFinder` (`offtarget_finder.py`)

- Implement genomic grouping (by `sseqid` / Chromosome).
- Implement a "Sliding Window" pairing logic:
  - For each Forward hit, find all Reverse hits on the same strand within $N$ bp.
- Classify amplicons into:

  - **On-Target**: Matches the intended junction coordinates.
  - **Cross-Target**: Formed by two primers from the same panel but different junctions.
  - **Off-Target**: Formed by one or more primers binding to an unintended genomic location.

### Phase 3: Integrate `ntthal`

- Utilize `primer3-py`'s `calc_heterodimer` (which wraps `ntthal`) to score every binding site found by BLAST against its local genomic sequence.
- Add an `off_target_dg` field to the `PrimerPair` object.

### Phase 4: Scoring and Selection

- Update the `MultiplexCostFunction` to incorporate the `off_target_dg` and count of potential off-target amplicons.

## 5. Pros and Cons of Alternatives

| Alternative | Pros | Cons |
| :--- | :--- | :--- |
| **BWA/Bowtie2 Mapping** | Extremely fast; handles millions of primers. | Hard to integrate thermodynamic 3' bonus; complex dependency. |
| **Dedicated isPcr Binary** | Fast and standard. | Adds another external dependency; hard to customize scoring. |
| **Plexus Hybrid (BLAST + ntthal)** | **High precision**; uses existing dependencies; customizable logic. | Slower than pure mapping for very large panels. |

## 6. Implementation Milestones

1. **Drafting**: Finalize `ISPCR_PLAN.md` (Current).
2. **Refactoring**: Clean up `offtarget_finder.py` to handle pairing.
3. **Integration**: Link `ntthal` scoring to the pairing logic.
4. **Validation**: Verify against known off-target fixtures in `tests/data`.
