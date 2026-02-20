# ================================================================================
# Oligonucleotide thermodynamic alignment calculations.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2025 Stefan Filges
#
# Based on the oligotm.c library from Primer3.
# Copyright (c) 1996-2007 Whitehead Institute for Biomedical Research
# ================================================================================

import math
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import primer3
from loguru import logger

from plexus.utils.utils import gc_content

# Kelvin to Celsius conversion factor
T_KELVIN = 273.15

# In cal/(K·mol)
GAS_CONSTANT = 1.987


@dataclass
class TmResult:
    """
    Result of Tm calculation.

        Tm: melting temperature
        bound: fraction bound
        ddG: gibb's energy
    """

    Tm: float
    bound: float
    ddG: float


class OligotmError(Exception):
    """Base exception for oligotm errors."""

    pass


class InvalidSequenceError(OligotmError):
    """Raised when sequence is invalid or too short."""

    pass


class InvalidConcentrationError(OligotmError):
    """Raised when concentration values are invalid."""

    pass


# SantaLucia (1998) parameters
# Entropy values (DS) in 0.1 cal/(K·mol)
DS_PARAMS = {
    "AA": 222,
    "AC": 224,
    "AG": 210,
    "AT": 204,
    "AN": 224,
    "CA": 227,
    "CC": 199,
    "CG": 272,
    "CT": 210,
    "CN": 272,
    "GA": 222,
    "GC": 244,
    "GG": 199,
    "GT": 224,
    "GN": 244,
    "TA": 213,
    "TC": 222,
    "TG": 227,
    "TT": 222,
    "TN": 227,
    "NA": 168,
    "NC": 210,
    "NG": 220,
    "NT": 215,
    "NN": 220,
}

# Enthalpy values (DH) in 100 cal/mol
DH_PARAMS = {
    "AA": 79,
    "AC": 84,
    "AG": 78,
    "AT": 72,
    "AN": 72,
    "CA": 85,
    "CC": 80,
    "CG": 106,
    "CT": 78,
    "CN": 78,
    "GA": 82,
    "GC": 98,
    "GG": 80,
    "GT": 84,
    "GN": 80,
    "TA": 72,
    "TC": 82,
    "TG": 85,
    "TT": 79,
    "TN": 72,
    "NA": 72,
    "NC": 80,
    "NG": 78,
    "NT": 72,
    "NN": 72,
}

# Delta G values in cal/mol
DG_PARAMS = {
    "AA": 1000,
    "AC": 1440,
    "AG": 1280,
    "AT": 880,
    "AN": 880,
    "CA": 1450,
    "CC": 1840,
    "CG": 2170,
    "CT": 1280,
    "CN": 1450,
    "GA": 1300,
    "GC": 2240,
    "GG": 1840,
    "GT": 1440,
    "GN": 1300,
    "TA": 580,
    "TC": 1300,
    "TG": 1450,
    "TT": 1000,
    "TN": 580,
    "NA": 580,
    "NC": 1300,
    "NG": 1280,
    "NT": 880,
    "NN": 580,
}


def calc_thermodynamics(seq: str) -> tuple:
    """
    Calculate enthalpy (dh) and entropy (ds) for sequence using
    SantaLucia (1998) parameters.

    Args:
        seq - Sequence to analyze
    """
    dh = 0
    ds = 0

    # Select parameter tables
    h_table = DH_PARAMS
    s_table = DS_PARAMS

    # Symmetry correction
    if symmetry(seq):
        ds += 14

    # Terminal AT penalty
    if seq[0] in "AT":
        ds += -41
        dh += -23
    elif seq[0] in "CG":
        ds += 28
        dh += -1

    if seq[-1] in "AT":
        ds += -41
        dh += -23
    elif seq[-1] in "CG":
        ds += 28
        dh += -1

    # Calculate nearest-neighbor contributions
    for i in range(len(seq) - 1):
        pair = seq[i : i + 2]
        dh += h_table.get(pair, h_table["NN"])
        ds += s_table.get(pair, s_table["NN"])

    return dh, ds


def oligotm(
    seq: str,
    dna_conc: float = 50.0,
    salt_conc: float = 50.0,
    divalent_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
    annealing_temp: float = 60.0,
) -> TmResult:
    """
    Calculate melting temperature using SantaLucia nearest-neighbor method.

    Args:
        seq: DNA sequence
        dna_conc: DNA concentration (nM)
        salt_conc: Monovalent salt concentration (mM)
        divalent_conc: Divalent cation concentration (mM)
        dntp_conc: dNTP concentration (mM)
        dmso_conc: DMSO concentration (%)
        dmso_fact: DMSO correction factor (default 0.6)
        formamide_conc: Formamide concentration (mol/L)
        annealing_temp: Annealing temperature for binding calculation

    Returns:
        TmResult with Tm and bound percentage

    Raises:
        InvalidSequenceError: If sequence is too short or contains invalid characters
        InvalidConcentrationError: If concentration values are invalid
    """
    if len(seq) < 2:
        raise InvalidSequenceError("Sequence must be at least 2 bases long")

    # Validate sequence contains only valid bases
    valid_bases = set("ATGCN")
    if not all(base in valid_bases for base in seq.upper()):
        raise InvalidSequenceError(
            "Sequence contains invalid characters. Only A, T, G, C, N allowed"
        )

    seq = seq.upper()

    # Validate and convert divalent to monovalent
    dv_to_mv = divalent_to_monovalent(divalent_conc, dntp_conc)

    # Count GC for formamide correction
    gc_count = sum(1 for base in seq if base in "GC")
    seq_len = len(seq)

    # Calculate thermodynamics
    dh, ds = calc_thermodynamics(seq)

    delta_H = dh * -100.0  # Convert to cal/mol
    delta_S = ds * -0.1  # Convert to cal/(K·mol)

    K_mM = salt_conc + dv_to_mv

    temp = annealing_temp + T_KELVIN

    # Apply salt correction (SantaLucia method)
    delta_S = delta_S + 0.368 * (seq_len - 1) * math.log(K_mM / 1000.0)

    # Calculate Gibbs free energy and equilibrium constant
    # Using salt-corrected delta_S as per primer3 C source
    ddG = delta_H - temp * delta_S
    ka = math.exp(-ddG / (GAS_CONSTANT * temp))

    if symmetry(seq):
        Tm = delta_H / (delta_S + GAS_CONSTANT * math.log(dna_conc / 1e9)) - T_KELVIN
        bound = (1 / (1 + math.sqrt(1 / ((dna_conc / 1e9) * ka)))) * 100
    else:
        Tm = delta_H / (delta_S + GAS_CONSTANT * math.log(dna_conc / 4e9)) - T_KELVIN
        bound = (1 / (1 + math.sqrt(1 / ((dna_conc / 4e9) * ka)))) * 100

    # Apply DMSO correction
    if dmso_conc > 0.0:
        Tm -= dmso_conc * dmso_fact

    # Apply formamide correction (independent of DMSO)
    if formamide_conc > 0.0:
        Tm += (0.453 * gc_count / seq_len - 2.88) * formamide_conc

    return TmResult(Tm=Tm, bound=bound, ddG=ddG)


def oligodg(seq: str) -> float:
    """
    Calculate delta G of oligonucleotide disruption using SantaLucia parameters.

    Raises:
        InvalidSequenceError: If sequence is too short
    """
    if len(seq) < 2:
        raise InvalidSequenceError("Sequence must be at least 2 bases long")

    seq = seq.upper()
    dg = -1960  # Initial dG

    # Terminal AT penalty
    if seq[0] in "AT":
        dg += -50
    if seq[-1] in "AT":
        dg += -50

    # Calculate nearest-neighbor contributions
    for i in range(len(seq) - 1):
        pair = seq[i : i + 2]
        dg += DG_PARAMS.get(pair, DG_PARAMS["NN"])

    if symmetry(seq):
        dg += -430  # Symmetry correction

    return dg / 1000.0


# Return the delta G of the last len bases of oligo if oligo is at least len
# bases long; otherwise return the delta G of oligo.
def end_oligodg(seq: str, length: int) -> float:
    """
    Calculate delta G for last 'length' bases of sequence. Return the delta G of the
    last len bases of oligo if oligo is at least len bases long; otherwise return
    the delta G of oligo.

    """
    seq_len = len(seq)
    if seq_len < length:
        return oligodg(seq)
    return oligodg(seq[-length:])


def symmetry(seq: str) -> bool:
    """Check if sequence is self-complementary/symmetrical."""
    seq_len = len(seq)
    if seq_len % 2 == 1:
        return False

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    mid = seq_len // 2

    for i in range(mid):
        if seq[i] not in complement or seq[-(i + 1)] != complement[seq[i]]:
            return False

    return True


def divalent_to_monovalent(divalent: float, dntp: float) -> float:
    """Convert divalent salt concentration to monovalent equivalent."""
    if divalent == 0:
        dntp = 0
    if divalent < 0 or dntp < 0:
        raise InvalidConcentrationError(
            "Divalent and dNTP concentrations must be non-negative"
        )
    if divalent < dntp:
        divalent = dntp
    return 120 * math.sqrt(divalent - dntp)


# Both functions return the melting temperature of the given oligo
# calculated as specified by user, but oligotm _should_ only be used on
# DNA sequences of length <= MAX_PRIMER_LENGTH (which is defined
# elsewhere).  seqtm uses oligotm for sequences of length <=
# MAX_PRIMER_LENGTH, and a different, G+C% based formula for longer
# sequences.  For oligotm(), no error is generated on sequences
# longer than MAX_PRIMER_LENGTH, but the formula becomes less
# accurate as the sequence grows longer.
def long_seq_tm(
    seq: str,
    start: int,
    length: int,
    salt_conc: float,
    divalent_conc: float,
    dntp_conc: float,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> TmResult:
    """
    Calculate Tm for long sequences using GC% formula.

    Uses: Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length

    Raises:
        InvalidSequenceError: If sequence bounds are invalid
        InvalidConcentrationError: If concentration values are invalid
    """
    dv_to_mv = divalent_to_monovalent(divalent_conc, dntp_conc)
    salt_conc = salt_conc + dv_to_mv

    if start < 0 or length <= 0 or start + length > len(seq):
        raise InvalidSequenceError(
            f"Invalid sequence bounds: start={start}, length={length}, seq_len={len(seq)}"
        )

    subseq = seq[start : start + length].upper()
    gc_count = sum(1 for base in subseq if base in "GC")

    Tm = (
        81.5
        - dmso_conc * dmso_fact
        + (0.453 * gc_count / length - 2.88) * formamide_conc
        + 16.6 * math.log10(salt_conc / 1000.0)
        + 41.0 * (gc_count / length)
        - 600.0 / length
    )

    # GC% formula has no thermodynamic ddG; use 0.0 (matches bound=0.0 convention)
    return TmResult(Tm=Tm, bound=0.0, ddG=0.0)


def seqtm(
    seq: str,
    dna_conc: float = 50.0,
    salt_conc: float = 50.0,
    divalent_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
    nn_max_len: int = 60,
    annealing_temp: float = 60.0,
) -> TmResult:
    """
    Calculate melting temperature for sequences of any length.

    Uses SantaLucia nearest-neighbor model for sequences <= nn_max_len,
    GC% formula for longer sequences.

    Args:
        nn_max_len - Max Value for which nearest neighbour model is valid. Default 60.
    """
    seq_len = len(seq)

    if seq_len > nn_max_len:
        return long_seq_tm(
            seq,
            0,
            seq_len,
            salt_conc,
            divalent_conc,
            dntp_conc,
            dmso_conc,
            dmso_fact,
            formamide_conc,
        )
    else:
        return oligotm(
            seq,
            dna_conc,
            salt_conc,
            divalent_conc,
            dntp_conc,
            dmso_conc,
            dmso_fact,
            formamide_conc,
            annealing_temp,
        )


def calculate_single_primer_thermodynamics(primer_list, config, orientation: str):
    """
    Function to calculate thermodynamic properties of a list of individual primers.
    Sequence-specific parameters (e.g. GC-content, GC clamp, etc.) are calculated
    elsewhere. Primers for this function should already have adequate sequence
    properties to avoid calculations on primers that would fail other criteria.

    Args:
        primer_list - List of primer sequences, e.g. LEFT_PRIMERS or RIGHT_PRIMERS.
        config - Configuration file containing decision thresholds.
    """
    # primer counter
    primers_considered = 0

    # Exclusion counters
    # Order of filtering is important as only the first reason to fail will be
    # incremented by 1, all subsequent tests won't be run. E.g. a primer may fail
    # both tm_too_low and self_dimer_any but tm_too_low is checked first.
    # tm_too_low will be raised by 1, self_dimer_any will not be checked and
    # stays the same.
    tm_too_low = 0
    tm_too_high = 0
    amount_bound_low = 0
    amount_bound_high = 0
    self_dimer_any = 0
    self_dimer_end = 0
    self_hairpin = 0
    primer_3prime_too_stable = 0

    # Get config sections for easier access
    singleplex = config.singleplex_design_parameters
    pcr = config.pcr_conditions

    max_tm = singleplex.PRIMER_MAX_TM
    min_tm = singleplex.PRIMER_MIN_TM

    max_bound = singleplex.PRIMER_MAX_BOUND
    min_bound = singleplex.PRIMER_MIN_BOUND

    primer_max_hairpin_tm = singleplex.PRIMER_MAX_HAIRPIN_TH
    primer_max_self_any = singleplex.PRIMER_MAX_SELF_ANY_TH
    primer_max_self_end = singleplex.PRIMER_MAX_SELF_END_TH
    primer_max_end_stability = singleplex.PRIMER_MAX_END_STABILITY

    good_primers = []

    # Iterate over all primers in the list. Hard filters are ordered such that more
    # common reasons for exclusion are first, so that fewer calculations have
    # to be made. Penalties are calculated only for primers that passed
    # all hard filters.
    for primer in primer_list:
        primers_considered += 1

        # Define ThermoAnalysis object and set parameters
        oligo_calc = primer3.thermoanalysis.ThermoAnalysis()

        oligo_calc.set_thermo_args(
            mv_conc=pcr.mv_concentration,
            dv_conc=pcr.dv_concentration,
            dntp_conc=pcr.dntp_concentration,
            dna_conc=pcr.primer_concentration,
            dmso_conc=pcr.dmso_concentration,
            dmso_fact=pcr.dmso_fact,
            formamide_conc=pcr.formamide_concentration,
            annealing_temp_c=pcr.annealing_temperature,
            temp_c=37.0,
            max_nn_length=60,
            max_loop=30,
            tm_method="santalucia",
            salt_corrections_method="santalucia",
        )

        # Initial primer penalty
        primer_penalty = 0

        # Sanity check primer size, these should have been handeled elsewhere already
        primer_length = len(primer.seq)
        if (primer_length > singleplex.primer_max_length) or (
            primer_length < singleplex.primer_min_length
        ):
            raise ValueError(
                f"Primer length of {primer_length} does not match thresholds."
            )

        # Sanity check: PRIMER_MAX_GC and PRIMER_MIN_GC
        primer_gc = gc_content(primer.seq)
        if (primer_gc > singleplex.primer_max_gc) or (
            primer_gc < singleplex.primer_min_gc
        ):
            raise ValueError(f"GC-content of {primer_gc} does not match thresholds.")

        # Primer Melting temperature and amount bound
        tm_bound = seqtm(
            seq=primer.seq,
            salt_conc=pcr.mv_concentration,
            divalent_conc=pcr.dv_concentration,
            dntp_conc=pcr.dntp_concentration,
            dna_conc=pcr.primer_concentration,
            dmso_conc=pcr.dmso_concentration,
            dmso_fact=pcr.dmso_fact,
            formamide_conc=pcr.formamide_concentration,
            annealing_temp=pcr.annealing_temperature,
        )
        primer_tm = round(tm_bound.Tm, 1)
        primer_bound = round(tm_bound.bound, 1)
        # Alternative to calculate Tm using primer3-py, should be identical to seqtm,
        # but does not produce amount bound
        # primer_tm = round(oligo_calc.calc_tm(primer),1)

        # Check: PRIMER_MIN_TM and PRIMER_MAX_TM
        if primer_tm > max_tm:
            tm_too_high += 1
            continue
        if primer_tm < min_tm:
            tm_too_low += 1
            continue

        # Check: PRIMER_MIN_BOUND and PRIMER_MAX_BOUND
        # Primer3 defaults are PRIMER_MAX_BOUND = 110 and PRIMER_MIN_BOUND = -10, which
        # would never fail this check!
        #
        # SantaLucia argues that primers should not be matched on melting temperature
        # (PRIMER_OPT_TM) but on the fraction of primers bound at annealing temperature
        # (PRIMER_OPT_BOUND). Especially multiplex primers should profit from
        # thermodynamic parameters where the individual primers match better
        # to each other.
        if primer_bound > max_bound:
            amount_bound_high += 1
            continue
        if primer_bound < min_bound:
            amount_bound_low += 1
            continue

        # PRIMER_MAX_HAIRPIN_TH
        # This is the most stable monomer structure of internal oligo calculated by
        # the thermodynamic approach.
        PRIMER_MAX_HAIRPIN_TH = oligo_calc.calc_hairpin(primer.seq).tm
        if PRIMER_MAX_HAIRPIN_TH > primer_max_hairpin_tm:
            self_hairpin += 1
            continue

        # Check: PRIMER_MAX_SELF_ANY_TH
        # The melting temperature of the most stable structure is calculated.
        PRIMER_SELF_ANY_TH = oligo_calc.calc_homodimer(primer.seq).tm
        if PRIMER_SELF_ANY_TH > primer_max_self_any:
            self_dimer_any += 1
            continue

        # PRIMER_MAX_SELF_END_TH
        # tries to bind the 3'-END to an identical primer and scores the best binding it
        # can find.
        PRIMER_MAX_SELF_END_TH = oligo_calc.calc_end_stability(
            primer.seq, primer.seq
        ).tm
        if PRIMER_MAX_SELF_END_TH > primer_max_self_end:
            self_dimer_end += 1
            continue

        # PRIMER_END_STABILITY
        # The maximum stability for the last five 3' bases of a left or right primer.
        # Bigger numbers mean more stable 3' ends.
        PRIMER_END_STABILITY = end_oligodg(primer.seq, 5)
        if PRIMER_END_STABILITY > primer_max_end_stability:
            primer_3prime_too_stable += 1
            continue

        # logger.info(f'{primers_considered} - Primer: {primer}')
        # logger.info(f'Primer Tm: {primer_tm}')
        # logger.info(f'Amount bound: {primer_bound}')
        # logger.info(f'PRIMER_MAX_SELF_ANY_TH: {round(PRIMER_SELF_ANY_TH,1)}')
        # logger.info(f'PRIMER_MAX_HAIRPIN_TH: {round(PRIMER_MAX_HAIRPIN_TH,1)}')
        # logger.info(f'PRIMER_MAX_SELF_END_TH: {round(PRIMER_MAX_SELF_END_TH,1)}')
        # logger.info(f'PRIMER_MAX_END_STABILITY: {round(PRIMER_END_STABILITY,1)}')

        # ==============================================================================
        # Now only "good" primers are left and we calculate the penalty for the primer
        # Similar to primer3: https://primer3.org/manual.html#calculatePenalties
        # ==============================================================================

        # Optimal values
        PRIMER_OPT_SIZE = singleplex.PRIMER_OPT_SIZE
        PRIMER_OPT_TM = (
            singleplex.PRIMER_OPT_TM
        )  # does not count if PRIMER_WT_TM_GT and PRIMER_WT_TM_LT are 0
        PRIMER_OPT_BOUND = (
            singleplex.PRIMER_OPT_BOUND
        )  # does not count if PRIMER_WT_BOUND_GT and PRIMER_WT_BOUND_LT are 0
        PRIMER_OPT_GC_PERCENT = (
            singleplex.PRIMER_OPT_GC_PERCENT
        )  # does not count if PRIMER_WT_GC_PERCENT_GT and PRIMER_WT_GC_PERCENT_LT are 0

        # Primer penalty weights
        PRIMER_WT_SIZE_LT = singleplex.PRIMER_WT_SIZE_LT
        PRIMER_WT_SIZE_GT = singleplex.PRIMER_WT_SIZE_GT
        PRIMER_WT_TM_GT = singleplex.PRIMER_WT_TM_GT
        PRIMER_WT_TM_LT = singleplex.PRIMER_WT_TM_LT
        PRIMER_WT_BOUND_GT = singleplex.PRIMER_WT_BOUND_GT
        PRIMER_WT_BOUND_LT = singleplex.PRIMER_WT_BOUND_LT
        PRIMER_WT_GC_PERCENT_GT = singleplex.PRIMER_WT_GC_PERCENT_GT
        PRIMER_WT_GC_PERCENT_LT = singleplex.PRIMER_WT_GC_PERCENT_LT
        PRIMER_WT_SELF_ANY_TH = singleplex.PRIMER_WT_SELF_ANY_TH
        PRIMER_WT_SELF_END_TH = singleplex.PRIMER_WT_SELF_END_TH
        PRIMER_WT_HAIRPIN_TH = singleplex.PRIMER_WT_HAIRPIN_TH
        PRIMER_WT_END_STABILITY = singleplex.PRIMER_WT_END_STABILITY

        tm_adjusted = primer_tm - 5

        if primer_tm > PRIMER_OPT_TM:
            primer_penalty += PRIMER_WT_TM_GT * (primer_tm - PRIMER_OPT_TM)
        if primer_tm < PRIMER_OPT_TM:
            primer_penalty += PRIMER_WT_TM_LT * (PRIMER_OPT_TM - primer_tm)
        if primer_bound > PRIMER_OPT_BOUND:
            primer_penalty += PRIMER_WT_BOUND_GT * (primer_bound - PRIMER_OPT_BOUND)
        if primer_bound < PRIMER_OPT_BOUND:
            primer_penalty += PRIMER_WT_BOUND_LT * (PRIMER_OPT_BOUND - primer_bound)
        if primer_gc > PRIMER_OPT_GC_PERCENT:
            primer_penalty += PRIMER_WT_GC_PERCENT_GT * (
                primer_gc - PRIMER_OPT_GC_PERCENT
            )
        if primer_gc < PRIMER_OPT_GC_PERCENT:
            primer_penalty += PRIMER_WT_GC_PERCENT_LT * (
                PRIMER_OPT_GC_PERCENT - primer_gc
            )
        if primer_length > PRIMER_OPT_SIZE:
            primer_penalty += PRIMER_WT_SIZE_GT * (primer_length - PRIMER_OPT_SIZE)
        if primer_length < PRIMER_OPT_SIZE:
            primer_penalty += PRIMER_WT_SIZE_LT * (PRIMER_OPT_SIZE - primer_length)

        # Self-complementarity penalty
        if tm_adjusted <= PRIMER_SELF_ANY_TH:
            primer_penalty += PRIMER_WT_SELF_ANY_TH * (
                PRIMER_SELF_ANY_TH - (tm_adjusted - 1)
            )
        else:
            primer_penalty += PRIMER_WT_SELF_ANY_TH * (
                1 / (tm_adjusted + 1 - PRIMER_SELF_ANY_TH)
            )

        # Self-end complementarity penalty
        if tm_adjusted <= PRIMER_MAX_SELF_END_TH:
            primer_penalty += PRIMER_WT_SELF_END_TH * (
                PRIMER_MAX_SELF_END_TH - (tm_adjusted - 1)
            )
        else:
            primer_penalty += PRIMER_WT_SELF_END_TH * (
                1 / (tm_adjusted + 1 - PRIMER_MAX_SELF_END_TH)
            )

        # Hairpin penalty
        if tm_adjusted <= PRIMER_MAX_HAIRPIN_TH:
            primer_penalty += PRIMER_WT_HAIRPIN_TH * (
                PRIMER_MAX_HAIRPIN_TH - (tm_adjusted - 1)
            )
        else:
            primer_penalty += PRIMER_WT_HAIRPIN_TH * (
                1 / (tm_adjusted + 1 - PRIMER_MAX_HAIRPIN_TH)
            )

        # Primer 3'-end stability penalty
        primer_penalty += PRIMER_WT_END_STABILITY * PRIMER_END_STABILITY

        # Add values to Primer object
        primer.bound = round(primer_bound, 2)
        primer.tm = round(primer_tm, 2)
        primer.gc = round(primer_gc, 2)
        primer.penalty = round(primer_penalty, 1)
        primer.self_any_th = round(PRIMER_SELF_ANY_TH, 2)
        primer.self_end_th = round(PRIMER_MAX_SELF_END_TH, 2)
        primer.hairpin_th = round(PRIMER_MAX_HAIRPIN_TH, 2)
        primer.end_stability = round(PRIMER_END_STABILITY, 2)
        primer.engine = "custom"

        # Save good primers
        good_primers.append(primer)

    # Text-based summary of the design process, only use non-zero criteria.
    eval_string = f"considered: {primers_considered}, "

    if tm_too_high > 0:
        eval_string += f"tm (high) > {max_tm}: {tm_too_high}, "
    if tm_too_low > 0:
        eval_string += f"tm (low) < {min_tm}: {tm_too_low}, "
    if amount_bound_low > 0:
        eval_string += f"Amount bound < {min_bound}: {amount_bound_low}, "
    if amount_bound_high > 0:
        eval_string += f"Amount bound > {max_bound}: {amount_bound_high}, "
    if self_hairpin > 0:
        eval_string += f"self hairpin tm > {primer_max_hairpin_tm}: {self_hairpin}, "
    if self_dimer_any > 0:
        eval_string += f"any self dimer tm > {primer_max_self_any}: {self_dimer_any}, "
    if self_dimer_end > 0:
        eval_string += f"end self dimer tm > {primer_max_self_end}: {self_dimer_end}, "
    if primer_3prime_too_stable > 0:
        eval_string += (
            f"3'-end stability > {primer_max_end_stability}: {primer_3prime_too_stable}"
        )

    logger.info(eval_string)

    return (good_primers, eval_string)


def calculate_single_primer_thermodynamics_parallel(left_kmers, right_kmers, config):
    """
    Multithreaded version of 'calculate_single_primer_thermodynamics'.
    """
    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor(
        max_workers=max(2, multiprocessing.cpu_count() - 2)
    ) as executor:
        # Submit both tasks
        left_future = executor.submit(
            calculate_single_primer_thermodynamics, left_kmers, config, "left"
        )
        right_future = executor.submit(
            calculate_single_primer_thermodynamics, right_kmers, config, "right"
        )

        return left_future.result(), right_future.result()
