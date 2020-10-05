# pdb-profiling

[![Build](https://img.shields.io/travis/naturegeorge/pdb-profiling?style=for-the-badge&logo=travis)](https://github.com/naturegeorge/pdb-profiling)
[![SupportPythonVersion](https://img.shields.io/pypi/pyversions/pdb-profiling.svg?style=for-the-badge&logo=python)](https://pypi.org/project/pdb-profiling/)
[![Version](https://img.shields.io/pypi/v/pdb-profiling?style=for-the-badge&logo=PYPI)](https://github.com/naturegeorge/pdb-profiling/blob/master/pdb_profiling/__init__.py)
[![Dependencies](https://img.shields.io/librariesio/github/NatureGeorge/pdb-profiling?style=for-the-badge&logo=PYPI)](https://github.com/naturegeorge/pdb-profiling/blob/master/setup.py)
[![PyPIDownloads](https://img.shields.io/pypi/dm/pdb-profiling.svg?style=for-the-badge&logo=PYPI)](https://pypi.org/project/pdb-profiling/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge&logo=github)](https://github.com/naturegeorge/pdb-profiling/blob/master/LICENSE)
[![GitHubDownloads](https://img.shields.io/github/downloads/NatureGeorge/pdb-profiling/total?style=for-the-badge&logo=github)](https://github.com/NatureGeorge/pdb-profiling/releases/)

![cover](https://user-images.githubusercontent.com/43134199/95018149-58cfc200-0690-11eb-9e64-760faec5130f.png)

Profiling Protein Structures from Protein Data Bank and integrate various resources.

[![HitCount](http://hits.dwyl.com/naturegeorge/pdb-profiling.svg)](http://hits.dwyl.com/naturegeorge/pdb-profiling)

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
  * UniProt API
    * <https://www.uniprot.org/uploadlists/>
    * <https://www.uniprot.org/uniprot/*.fasta>
  * Interactome3D API
    * <https://interactome3d.irbbarcelona.org/>
  * ModBase API (?)
  * Ensembl REST API
    * <https://rest.ensembl.org/documentation>
    * <https://rest.ensembl.org/documentation/info/sequence_id>
    * <https://rest.ensembl.org/documentation/info/archive_id_get>
  * Eutils API
    * <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/>
    * NOTE: currently only support minimum use
* Download data from PDB Archive against unexpected needs
  * wwPDB&RCSB: <https://ftp.wwpdb.org/pub/pdb/data/structures/>
  * EBI: <http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/>
  * wwPDB Versioned: <http://ftp-versioned.wwpdb.org/pdb_versioned/data/entries/>

## Install

```bash
pip install pdb-profiling
```

## Examples

See `examples/...`

1. [Introduction](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/Introduction.ipynb)
2. [Batch](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/Batch.ipynb)
3. [DisplayPDB](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/DisplayPDB.ipynb)
3. ...

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).