# Multiplex designer user guide


# Object classes

## Junction

Contains information for a specific junction.


## Primer

Contains information for a single primer.


## PrimerPair


## MultiplexPanel

The primary object.



## Problem 

We designed N primer pairs for each of k targets, meaning there are
$2*k*N$ available single primers. We now want to pick exactly one primer pair
flanking each of the k targets.

For each primer (of a primer pair) we need to calculate the interactions between
itself and any of the other primers in a potential multiplex. 

How do we choose the most optimal set of primers?


## On primer binding

SantaLucia argues, that primers should not be matched on melting temperature (PRIMER_OPT_TM) but on the fraction of primers bound at annealing temperature (PRIMER_OPT_BOUND). Especially multiplex primers should profit from thermodynamic parameters where the individual primers match better to each other. See also the chapter on [primer binding in the primer3 documentation](https://www.primer3plus.com/primer3plusHelp.html#primerBinding).

Below are two examples for configurations based on either Tm or fraction bound. 

**NOTE:** Fraction bound tend to lead to significantly fewer designs.

### Primer design based on Tm (classic)

Using the configuration below, primer are selected based on melting temperature (Tm). The permissible 
ranges for amount bound are such that primers are guaranteed to pass this criterion. Penalty weights
for amount bound are set to zero, i.e. deviation from the optimal amount bound does not receive a
penalty.

```
"PRIMER_OPT_TM": 60.0
"PRIMER_MIN_TM": 57.0
"PRIMER_MAX_TM": 63.0
"primer_pair_max_tm_difference": 3.0
"PRIMER_OPT_BOUND": 97.0
"PRIMER_MIN_BOUND": -10.0
"PRIMER_MAX_BOUND": 110.0
"PRIMER_WT_TM_GT": 1.0,
"PRIMER_WT_TM_LT": 1.0,
"PRIMER_WT_BOUND_GT": 0.0,
"PRIMER_WT_BOUND_LT": 0.0,
```

### Primer design based on fraction bound

This is eseentialy the opposite of the configuration for Tm-based-design. Deviation from optimal Tm
is not penalized, and permissible Tm ranges are such that most primers will pass these. Instead
fraction bound is used for selection.

```
"PRIMER_OPT_TM": 60.0
"PRIMER_MIN_TM": 40.0
"PRIMER_MAX_TM": 80.0
"primer_pair_max_tm_difference": 10.0
"PRIMER_OPT_BOUND": 97.0
"PRIMER_MIN_BOUND": 90.0
"PRIMER_MAX_BOUND": 110.0
"PRIMER_WT_TM_GT": 0.0,
"PRIMER_WT_TM_LT": 0.0,
"PRIMER_WT_BOUND_GT": 1.0,
"PRIMER_WT_BOUND_LT": 1.0,
```


