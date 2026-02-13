# Implementation specs for plexus

A package to design multiplex PCR primers.

## Overview

The user provides a list of target junctions (positions in the genome) and the package designs primers around these junctions. Sensible defaults are used for the parameters embedded in the package config, but can be altered by the user by providing a config file.

Junction file + config file -> design candidate primers -> check for SNPs -> identify cross dimers through thermodynamic alignment -> blast to identify off-targets -> optimize the multiplex -> output multiplex panel

## Modules

### designer

The designer module design the candidate primers for the given junctions.
