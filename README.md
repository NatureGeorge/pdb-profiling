# pdb-profiling

[![Version](https://img.shields.io/pypi/v/pdb-profiling?style=for-the-badge)](https://github.com/naturegeorge/pdb-profiling/blob/master/pdb_profiling/__init__.py)
[![Downloads](https://img.shields.io/pypi/dm/pdb-profiling?style=for-the-badge)](https://pypi.org/project/pdb-profiling/)
[![Dependences](https://img.shields.io/david/NatureGeorge/pdb-profiling?style=for-the-badge)](https://github.com/naturegeorge/pdb-profiling/blob/master/setup.py)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge)](https://github.com/naturegeorge/pdb-profiling/blob/master/LICENSE)

Profiling Protein Structures from Protein Data Bank and integrate various resources.

## Goal

* Gather helpful/insightful indexes to evaluate a PDB structure's usefulness in:
  * Entry level
  * Assembly level
  * Model level
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
  * PDBe REST API
    * <https://www.ebi.ac.uk/pdbe/api/doc/pdb.html>
    * <https://www.ebi.ac.uk/pdbe/api/doc/pisa.html>
    * <https://www.ebi.ac.uk/pdbe/api/doc/sifts.html>
  * PDBe Graph API (Neo4j Graph DataBase)
    * <https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/>
  * PDBe CoordinateServer API
    * <https://www.ebi.ac.uk/pdbe/coordinates/index.html>
  * PDBe ModelServer API
    * <https://www.ebi.ac.uk/pdbe/model-server/>
  * SWISS-MODEL Repository API
    * <https://swissmodel.expasy.org/docs/smr_openapi>
  * EBI Proteins API
    * <https://www.ebi.ac.uk/proteins/api/doc/>
  * Interactome3D API
    * <https://interactome3d.irbbarcelona.org/>
  * ModBase API (?)
* Release this project as a Python package

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).