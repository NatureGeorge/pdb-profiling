# How to define a representative set of protein structures?

## Goal

* Reduce redundant structures (Minimize redundancy)
* Cover most part of the protein (Maximize coverage)
* Select the structures with good quality
* Select those structures with biological significance and worth studying
* Multi criteria sorting and filtering
* Combine empirical threshold and data based threshold
* Refine to residue level

Let's first take a look at what predecessors have achieved.

## Enlarge representative set of protein structures

> UWE HOBOHM AND CHRIS SANDER _European Molecular Biology Laboratory, 69012 Heidelberg, Germany_ (1993)

### Criteria/Quality Control

* an initial filter
* a "softexclude" flag for marginal data
* a quality index $Q$
  * $Q=r+R/20$
  * $r$: resolution
  * $R$: R-factor (percent)
  * the one with higher Q is considered to be of lower quality.

#### Initial Filter

eliminate all chains with:

* `100%` sequence identity on the entire length to another chain, but `lower quality`
* more than `5%` of `nonstandard` or `unknown amino acids(UNK)`
* a length of less than `30` residues
* a resolution worse than `3.5`$A^{\circ}$
* an R-factor of more than `30%`
* models not based directly on X-ray or NMR data

> The values for cutoffs were chosen by experienceand can be changed
to meet special requirements.

#### Softexclude Flag

> Some chains are flagged for preferred exclusion inthe
subsequent selection procedure (“softexclude flag”). Such chains
remain in the list only if no homologous chain exists

These are chains for which:

* the number of residues with __sidechain coordinates__ is less than `90%` of the sequence length (e.g., backbone-only structures)
* the number of residues with __backbone coordinates__ is less than `90%` of the sequence length (e.g., C-α-only structures)
* the number of __alanine__ plus __glycine__ residues is higher than `40%` of the sequence length (e.g., structures with unknown sequence modeled as polyalanine)
* no data for resolution or R-factor are available

### Algorithms/Selection procedure

#### Algorithm 2 of Hobohm et al. (1992)

##### Feature

* removes redunant protein chains one by one
* follow "greedy" strategy

##### Procedure

* the chain with the largest number of neighbors is removed, until on neighbors are left
  * neighbors are here defined as pairs of chains with sequence identity above 25% (or for alignments shorter than `80` residues, with a sequence identity above the significance threshold derived by Sander & Schneider [1991])
* 
