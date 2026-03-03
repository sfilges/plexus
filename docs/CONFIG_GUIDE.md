# Configuration Guide for Plexus

Plexus uses a hierarchical JSON configuration to control every aspect of the multiplex PCR primer design process. This guide explains the parameters, their impact on the design, and provides sensible ranges for typical applications.

## Configuration Structure

The configuration is divided into six main sections:

1. `singleplex_design_parameters`: Individual primer properties.
2. `primer_pair_parameters`: Properties of the forward/reverse pair (e.g., amplicon size).
3. `pcr_conditions`: Thermodynamic environment (salt, concentrations).
4. `snp_check_parameters`: How variants (SNPs) affect primer scoring.
5. `blast_parameters`: Specificity and off-target detection thresholds.
6. `multiplex_picker_parameters`: Weights for the final panel optimization.

---

## 1. Singleplex Design Parameters

These parameters define what makes a "good" individual primer.

| Parameter | Default | Sensible Range | Description |
| :--- | :--- | :--- | :--- |
| `PRIMER_OPT_TM` | 60.0 | 58.0 - 62.0 | The ideal melting temperature (°C). |
| `PRIMER_MIN_TM` | 57.0 | 55.0 - 59.0 | Minimum allowed Tm. |
| `PRIMER_MAX_TM` | 63.0 | 61.0 - 68.0 | Maximum allowed Tm. |
| `PRIMER_OPT_SIZE` | 22 | 18 - 26 | The ideal primer length (bp). |
| `primer_min_length`| 18 | 15 - 20 | Minimum primer length. |
| `primer_max_length`| 28 | 25 - 35 | Maximum primer length. |
| `primer_min_gc` | 30 | 20 - 40 | Minimum GC content (%). |
| `primer_max_gc` | 70 | 60 - 80 | Maximum GC content (%). |
| `primer_gc_clamp` | 1 | 0 or 1 | If 1, requires 1-3 G/C bases in the last 5 bases of the 3' end. |
| `primer_max_poly_x`| 5 | 3 - 6 | Max allowed homopolymer length (e.g., AAAAA). |
| `PRIMER_MAX_HAIRPIN_TH` | 24.0 | 20.0 - 30.0 | Max allowed Tm (°C) for internal hairpin structures. |
| `PRIMER_MAX_SELF_ANY_TH` | 45.0 | 35.0 - 50.0 | Max allowed Tm (°C) for any self-dimer. |
| `PRIMER_MAX_SELF_END_TH` | 35.0 | 25.0 - 40.0 | Max allowed Tm (°C) for 3'-anchored self-dimers. |

---

## 2. Primer Pair Parameters

These parameters control the relationship between the forward and reverse primers.

| Parameter | Default | Sensible Range | Description |
| :--- | :--- | :--- | :--- |
| `PRIMER_PAIR_MAX_DIFF_TM` | 3.0 | 1.0 - 5.0 | Max difference in Tm between F and R primers. |
| `PRIMER_PRODUCT_OPT_SIZE` | 60 | 50 - 150 | Ideal total amplicon length (including primers). |
| `PRIMER_PRODUCT_MAX_SIZE` | 120 | 80 - 250 | Maximum total amplicon length. |
| `PRIMER_PRODUCT_MIN_INSERT_SIZE` | 20 | 10 - 50 | Minimum bases between the two primers. |

---

## 3. PCR Conditions

These values are used for thermodynamic calculations (SantaLucia 1998). Ensure these match your actual lab protocol for accurate Tm and Dimer predictions.

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `annealing_temperature` | 60.0 | The temperature used for the annealing step (°C). |
| `mv_concentration` | 50.0 | Monovalent cation concentration (mM), usually KCl. |
| `dv_concentration` | 1.5 | Divalent cation concentration (mM), usually MgCl2. |
| `dntp_concentration` | 0.6 | Total dNTP concentration (mM). |
| `primer_concentration` | 50.0 | Concentration of each individual primer (nM). |

---

## 4. SNP Check Parameters

Plexus queries VCF files to avoid designing primers over common variants.

| Parameter | Default | Sensible Range | Description |
| :--- | :--- | :--- | :--- |
| `af_threshold` | 0.01 | 0.001 - 0.05 | Minimum Allele Frequency (AF) to consider a SNP problematic. |
| `snp_penalty_weight` | 10.0 | 5.0 - 20.0 | Base penalty added to a primer for each overlapping SNP. |
| `snp_3prime_window` | 5 | 3 - 8 | Bases from the 3' end where a SNP is considered high-impact. |
| `snp_3prime_multiplier` | 3.0 | 2.0 - 5.0 | Multiplier for the penalty if the SNP is in the 3' window. |
| `snp_strict` | false | true/false | If true, any primer with a SNP > `af_threshold` is discarded immediately. |
| `snp_af_weight` | 0.0 | 0.0 - 1.0 | Exponent for AF scaling. `0.5` (sqrt) is recommended to penalize common SNPs more than rare ones. |

---

## 5. BLAST Parameters

Used to identify off-target products that might cause non-specific amplification.

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `length_threshold` | 15 | Min 3'-anchored alignment length (bp) to predict binding. |
| `max_mismatches` | 2 | Max mismatches allowed in a 3'-anchored alignment. |
| `max_amplicon_size` | 2000 | Max distance between two hits to be considered a potential off-target amplicon. |
| `ontarget_tolerance` | 5 | BP tolerance when verifying if a BLAST hit matches the intended target. |

---

## 6. Multiplex Picker Parameters

These weights control the final "cost" of a panel. The optimizer tries to minimize this cost.

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `wt_pair_penalty` | 1.0 | Weight for the basic primer quality penalty (Tm deviation, length, etc). |
| `wt_off_target` | 5.0 | Penalty for each predicted off-target product. Usually high to avoid mispriming. |
| `wt_cross_dimer` | 1.0 | Weight for penalties arising from primers in DIFFERENT pairs forming dimers. |
| `wt_pair_dimer` | 1.0 | Weight for the dimer score of the F/R primers within the SAME pair. |
| `wt_snp_penalty` | 3.0 | Weight for the accumulated SNP penalty from the `snpcheck` step. |
| `initial_solutions`| 100 | Number of iterations for stochastic selectors (Greedy, Simulated Annealing). |

---

## Best Practices

### For Small, High-Quality Panels (< 24 targets)

* Use `preset: "default"`.
* Set `snp_strict: true` to ensure no variants interfere with your clinical targets.
* Increase `initial_solutions` to `1000` for a more thorough search.

### For Large Discovery Panels (> 100 targets)

* Use `preset: "lenient"` to allow for a wider range of Tms.
* Keep `snp_strict: false` but use `snp_af_weight: 0.5` to prioritize the best available sites.
* Use the `Greedy` or `SimulatedAnnealing` selectors; avoid `BruteForce` or `DFS`.
