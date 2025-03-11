#-------------------------------------------// (3) OPTIMISER //----------------------------------
#
# ThermoBLAST is the algorithm that is used to scan groups of primer candidates against the human
# genome to deduce background false amplicons penalties (the penalties are a function of the percent 
# bound for each of the primers and the length of the false amplicon). ThermoBLAST is also used to 
# exclude primer candidates that cross-hybridize with desired amplicons. These two ThermoBLAST runs 
# are shown in solid green boxes in the flowchart figure above. The ThermoBLAST step maximizes 
# specificity by penalizing primers that have strong off-target hybridizations. See the white paper about 
# ThermoBLAST to learn more about how it works (download from www.dnasoftware.com).
# 
# Lastly, MultiPick performs a mixing and matching of FP and RP candidates from different targets 
# to make “multiplex solution sets”. MultiPick is run in 2 phases: first MultiPick analyzes all 100 
# solutions from all the targets to determine all the exclusions and to then form a shorter list of 
# 4 primer pair candidates for each of the N-plex targets. This shorter list of primer candidates 
# is submitted to ThermoBLAST against the human genome to deduce all the possible false amplicons. 
# For a 20-plex example, ThermoBLAST would consider 160 primers (i.e. 4 FP and 4 RP for 20 targets) 
# against the human genome – such massive runs are only possible due to the high-performance 
# cloud computing of Amazon Web Services and also clever algorithms written by DNA Software. 
# In phase 2, MultiPick considers all 4N possible multiplex combinations that are possible and 
# finds the combination(s) with the best score(s) consisting of Forward Primer and Reverse Primer 
# (and optionally Probe and/or Reverse Transcription primer). For large N-plexes, the exponential 
# explosion requires a sophisticated algorithm to solve the problem. DNA Software has developed s
# uch an algorithm – that uses a depth-first search with pruning algorithm, an approach used in
# 21st century artificial intelligence applications.
#-----------------------------------------------------------------------------------------------
