# ktools

ktools implements two varieties of Variational Transition State Theory (VTST):
Canonical (i.e. thermal) and J-resolved Microcanonical (i.e. for fixed internal energies and
angular momentum). From molecular vibration frequencies, moments of inertia, and potential
energy provided at points along a reaction path to compute variationally optimized sums of states
and J-resolved microcanonical rate constants (i.e. k(E,J)) for the reaction. These quantities can be
used in 2-dimensional (i.e. depending on both E and J) master equations. By summing over J, it
also gives microcanonical VTST rate constants (i.e. k(E) that can be used in 1-D (i.e. depending
on E, alone) master equations. An important feature otf the code is that it automatically identifies
cases where two or more bottlenecks occur along the same reaction path and then uses W. H.
Miller's unified statistical theory to compute the over-all effective rate constant.

ktools is part of the Multiwell program suite available here:

http://clasp-research.engin.umich.edu/multiwell/downloads.php

Current make file is configured to be placed in <multiwell-2017>/src/ktools and the binary is copied to <multiwell-2017>/bin
