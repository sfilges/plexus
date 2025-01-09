import primer3

# ---------------// Hairpin analysis //------------------
# Use primer3 to calculate hairpin
# max seq length is 60 nucleotides
# uses the thal library from primer3
# results stored as ThermoResult object

def calculate_hairpin(sequence, pcr_conditions):
    """
    Calculate hairpin thermodynamics for a DNA sequence.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_hairpin(
        seq=str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,  # Standard temperature for dG calculations
        max_loop=30,
        output_structure=False
    )

def calculate_melting_temperature(sequence, pcr_conditions):
    """
    Calculate melting temperature for a DNA sequence.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including concentrations and temperature
        
    Returns:
        float: Melting temperature in Celsius
    """
    return primer3.bindings.calc_tm(
        seq=str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        dmso_conc=pcr_conditions['dmso_concentration'],
        dmso_fact=pcr_conditions['dmso_fact'],
        formamide_conc=pcr_conditions['formamide_concentration'],
        annealing_temp_c=pcr_conditions['annealing_temperature'],
        max_nn_length=60,
        tm_method='santalucia',
        salt_corrections_method='santalucia'
    )

def calculate_homodimer(sequence, pcr_conditions):
    """
    Calculate homodimer formation thermodynamics.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_homodimer(
        str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,
        max_loop=30,
        output_structure=False
    )

def calculate_heterodimer(sequence1, sequence2, pcr_conditions):
    """
    Calculate heterodimer formation thermodynamics between two sequences.
    
    Args:
        sequence1 (str or Seq): First DNA sequence
        sequence2 (str or Seq): Second DNA sequence
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_heterodimer(
        str(sequence1),
        str(sequence2),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,
        max_loop=30,
        output_structure=False
    )