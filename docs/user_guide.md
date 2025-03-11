# Multiplex designer user guide




## Problem 

We designed N primer pairs for each of k targets, meaning there are
$2*k*N$ available single primers. We now want to pick exactly one primer pair
flanking each of the k targets.

For each primer (of a primer pair) we need to calculate the interactions between
itself and any of the other primers in a potential multiplex. 

How do we choose the most optimal set of primers?