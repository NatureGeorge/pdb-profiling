# pdb-profiling

Profiling Protein Structures from Protein Data Bank and integrate various resources.

## Goal

* Gather helpful/insightful indexes to evaluate a PDB structure's usefulness in:
  * Entry level
  * Model level (?)
  * Entity level
  * Chain level
  * Residue level
* Define the representative set of protein structures:
  * of a cluster with nearly identical sequences
  * of UniProt Entry
  * of UniProt Isoform
  * or any other assigned structure dataset
* Provide interface for ID/residue mapping
* Apply mature and robust API to collect well-organized data
  * Implement PDBe REST API (i.e <https://www.ebi.ac.uk/pdbe/api/doc/pdb.html>)
  * Implement PDBe Graph Database API (<https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/>)
* Release this project as a Python package

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).